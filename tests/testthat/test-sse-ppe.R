skip_on_cran()

test_that("the MLE recovers a known noncentrality parameter", {
  for (spec in list(list(df = 1, ncp = 8), list(df = 3, ncp = 15), list(df = 5, ncp = 25))) {
    d <- ppe_dofv(n = 4000L, df = spec$df, ncp = spec$ncp, seed = 7L)
    fit <- .ppeChiSquareMle(d, df = spec$df)
    expect_equal(fit$estimate, spec$ncp, tolerance = 0.05)
  }
})

test_that("the MLE recovers a known df when ncp is held at zero", {
  d <- ppe_dofv(n = 4000L, df = 4, ncp = 0, seed = 11L)
  fit <- .ppeChiSquareMle(d, ncp = 0)
  expect_equal(fit$estimate, 4, tolerance = 0.05)
})

test_that("the MLE agrees with a direct likelihood grid", {
  d <- ppe_dofv(n = 300L, df = 2, ncp = 6, seed = 13L)
  grid <- seq(0.1, 20, by = 0.01)
  ll <- vapply(grid, function(ncp) sum(stats::dchisq(d[d > 0], 2, ncp, log = TRUE)), numeric(1))

  expect_equal(.ppeChiSquareMle(d, df = 2)$estimate, grid[which.max(ll)], tolerance = 0.02)
})

test_that("non-positive values are counted, never hidden", {
  d <- c(-2, -1, 0, 3, 4, 5, 6)
  fit <- .ppeChiSquareMle(d, df = 1)

  expect_equal(fit$nRetained, 4L)
  expect_equal(fit$nNonPositive, 3L)
  expect_equal(fit$n, 7L)
})

test_that("a boundary solution is flagged rather than reported as a clean fit", {
  d <- c(0.1, 0.2, 0.15, 0.3)
  fit <- .ppeChiSquareMle(d, df = 4)

  expect_true(fit$boundary)
  expect_lt(fit$estimate, 1e-6)
})

test_that("the boundary warning's mode-specific detail line never claims interval degeneracy", {
  # Unit-level, deterministic coverage of the mode split: a df-target
  # boundary solution is numerically implausible to reach through the real
  # optimizer with continuous chi-square data (verified empirically -- even
  # 5000 retained values at .Machine$double.xmin converge to a df estimate
  # around 0.003, well above the 1e-8 boundary threshold), so this checks
  # the message-selection logic directly rather than searching for a seed
  # that may not exist.
  ncp_detail <- .ppeBoundaryDetail("ncp")
  df_detail <- .ppeBoundaryDetail("df")

  expect_true(grepl("power", ncp_detail, ignore.case = TRUE))
  expect_false(grepl("degenerate", ncp_detail, fixed = TRUE))
  expect_false(grepl("power", df_detail, ignore.case = TRUE))
  expect_false(grepl("degenerate", df_detail, fixed = TRUE))
})

test_that("ppeSummary's boundary warning does not claim the interval degenerates (power mode, end-to-end)", {
  # Power mode: ncp pinned at its lower bound. Seed 75 with df = 4, ncp = 0,
  # n = 60 reliably lands on boundary = TRUE (verified by direct search).
  sse_power <- fake_ppe_sse_object(df = 4, ncp = 0, n = 60L, seed = 75L)
  cmp_power <- sseComparison("simulation", "alt1", df = 4)
  msg_power <- tryCatch(
    { ppeSummary(sse_power, comparisons = cmp_power, bootstrapSamples = 20L, bootSeed = 1L); NA_character_ },
    warning = function(w) conditionMessage(w)
  )
  expect_false(is.na(msg_power))
  # The old false claim was "the interval is degenerate", stated as fact.
  # The corrected message may still mention degeneracy, but only hedged
  # ("does not generally degenerate") -- never asserted outright.
  expect_false(grepl("is degenerate", msg_power, fixed = TRUE))
  expect_true(grepl("does not generally degenerate", msg_power, fixed = TRUE))
  expect_true(grepl("power", msg_power, ignore.case = TRUE))
  expect_true(grepl("nonregular", msg_power, ignore.case = TRUE))
})

test_that("a boundary bootstrap interval actually has a positive upper endpoint, not just non-degenerate wording", {
  # Behavioral proof, not just a check on the warning's prose: at this
  # boundary solution (ncp pinned at its lower bound), finite bootstrap
  # samples drawn under ncp = 0 (df = 4) often refit to a strictly positive
  # noncentrality -- a direct 1000-sample reproduction found positive
  # refits in 45.5% of draws (see the technical reference's PPE section).
  # 300 bootstrap samples makes at least one positive refit overwhelmingly
  # likely (P(none positive) = (1 - 0.455)^300, astronomically small),
  # while still running quickly.
  sse_power <- fake_ppe_sse_object(df = 4, ncp = 0, n = 60L, seed = 75L)
  cmp_power <- sseComparison("simulation", "alt1", df = 4)

  s <- suppressWarnings(
    ppeSummary(sse_power, comparisons = cmp_power, bootstrapSamples = 300L, bootSeed = 1L)
  )

  expect_true(s$boundary)
  expect_equal(s$estimate, 0, tolerance = 1e-6)
  expect_gt(s$ci_upper, 0)
  expect_gt(s$ci_upper, s$ci_lower)
})

test_that("too few positive values is an error, not a silent fit", {
  expect_error(.ppeChiSquareMle(c(-1, 2), df = 1), "at least 2")
})

test_that(".ppeChiSquareMle rejects both df and ncp supplied together", {
  expect_error(.ppeChiSquareMle(c(1, 2, 3), df = 1, ncp = 2), "exactly one")
})

test_that(".ppeChiSquareMle rejects neither df nor ncp supplied", {
  # c(-1, 1) has only 1 positive value, so it discriminates the validation
  # order: with df/ncp checked first (as implemented), this reports "exactly
  # one" -- with the "at least 2 positive" check first, it would instead
  # report the (also true, but not the actual) problem: too few positive
  # values. c(1, 2, 3) does NOT discriminate this -- both orderings emit
  # "exactly one" for it, since 3 positive values pass the count check
  # regardless of order.
  expect_error(.ppeChiSquareMle(c(-1, 1)), "exactly one")
})

test_that("the nonpositive policy controls warning and abort, never the counts", {
  expect_warning(.ppeApplyNonpositivePolicy(2L, 5L, "warn"), "2 of 5")
  expect_error(.ppeApplyNonpositivePolicy(2L, 5L, "error"), "2 of 5")
  expect_silent(.ppeApplyNonpositivePolicy(2L, 5L, "drop"))
})

test_that("the nonpositive-policy warning does not claim a knowable bias direction", {
  # Regression test for the false "biased upward" truncation-renormalization
  # claim: P(X > 0) = 1 for a noncentral chi-square, so there is no missing
  # normalization constant, and the exclusion's bias direction is not
  # knowable from that fact alone.
  msg <- tryCatch(
    { .ppeApplyNonpositivePolicy(2L, 5L, "warn"); NA_character_ },
    warning = function(w) conditionMessage(w)
  )
  expect_false(grepl("biased upward", msg, fixed = TRUE))
  expect_true(grepl("selected-subset|selection", msg))
})

# --- .ppePowerPlotData(method = "distribution_mle") ---------------------

test_that("the fitted noncentrality does not depend on the threshold, unlike exceedance", {
  # This is the defining difference from method = "exceedance": one
  # noncentrality fitted from the whole distribution, reused for every
  # threshold, rather than a fresh implied noncentrality per threshold.
  #
  # df = 2 is passed explicitly via sseComparison(): fake_ppe_sse_object()'s
  # fixture model spec (fake_sse_fit_specs()) always implies df = 1 from its
  # own parameter counts, independent of the `df` used to generate the draws
  # -- so recovering the true generating ncp needs the true df supplied too.
  sse <- fake_ppe_sse_object(df = 2, ncp = 10, n = 500L, seed = 21L)
  cmp <- sseComparison("simulation", "alt1", df = 2)
  data <- .ppePowerPlotData(
    sse,
    thresholds = c(4, 12),
    studySizes = 12L,
    comparisons = cmp,
    nonpositivePolicy = "drop"
  )

  ncp_by_threshold <- split(data$mle_ncp, data$threshold)
  expect_length(unique(vapply(ncp_by_threshold, unique, numeric(1))), 1L)

  # And, unlike exceedance, that single fitted value should recover the
  # generating ncp (loosely -- this is a sanity check, not a tolerance test;
  # the tight tolerance checks already live in the .ppeChiSquareMle() tests).
  expect_equal(unique(data$mle_ncp), 10, tolerance = 0.1)
})

test_that("distribution_mle refuses a comparison with an explicit criticalValue", {
  sse <- fake_ppe_sse_object(df = 1, ncp = 8, n = 20L, seed = 5L)
  cmp <- sseComparison("simulation", "alt1", criticalValue = 5)

  expect_error(
    .ppePowerPlotData(sse, comparisons = cmp, nonpositivePolicy = "drop"),
    "not eligible"
  )
})

test_that("distribution_mle warns when df is inferred rather than declared", {
  # distribution_mle treats df as known -- it is the exponent in the
  # noncentral chi-square density the whole fit rests on -- so, unlike the
  # legacy exceedance branch, an inferred df must not be silent.
  sse <- fake_ppe_sse_object(df = 1, ncp = 8, n = 20L, seed = 5L)

  expect_warning(
    .ppePowerPlotData(sse, nonpositivePolicy = "drop"),
    "Degrees of freedom inferred"
  )
})

test_that("df_source reports explicit vs. inferred provenance", {
  sse <- fake_ppe_sse_object(df = 1, ncp = 8, n = 20L, seed = 5L)

  explicit_cmp <- sseComparison("simulation", "alt1", df = 1)
  explicit_data <- .ppePowerPlotData(
    sse,
    comparisons = explicit_cmp,
    nonpositivePolicy = "drop"
  )
  expect_true(all(explicit_data$df_source == "explicit"))

  inferred_data <- suppressWarnings(
    .ppePowerPlotData(sse, nonpositivePolicy = "drop")
  )
  expect_true(all(inferred_data$df_source == "parameter_count"))
})

# --- .ppeParametricBootstrap() -------------------------------------------

test_that("the bootstrap is reproducible under a fixed seed", {
  a <- .ppeParametricBootstrap(8, df = 1, nRetained = 100L,
                               bootstrapSamples = 50L, seed = 5L)
  b <- .ppeParametricBootstrap(8, df = 1, nRetained = 100L,
                               bootstrapSamples = 50L, seed = 5L)

  expect_equal(a$ci_lower, b$ci_lower)
  expect_equal(a$ci_upper, b$ci_upper)
})

test_that("the bootstrap interval brackets the point estimate", {
  res <- .ppeParametricBootstrap(8, df = 1, nRetained = 200L,
                                 bootstrapSamples = 200L, seed = 3L)

  expect_lt(res$ci_lower, 8)
  expect_gt(res$ci_upper, 8)
  expect_equal(res$n_successful, 200L)
  expect_equal(res$interval_type, "model_based")
})

test_that("the caller's RNG stream is left untouched", {
  set.seed(99L)
  before <- .Random.seed
  .ppeParametricBootstrap(8, df = 1, nRetained = 50L,
                          bootstrapSamples = 20L, seed = 5L)

  expect_equal(.Random.seed, before)
})

test_that("bootstrapSamples = 0 skips the bootstrap without changing the estimate", {
  res <- .ppeParametricBootstrap(8, df = 1, nRetained = 50L,
                                 bootstrapSamples = 0L, seed = 5L)

  expect_true(is.na(res$ci_lower))
  expect_equal(res$n_successful, 0L)
})

test_that("bootstrap failures are counted rather than aborting the fit", {
  # estimate = 1e-16, df = 1, nRetained = 3L (as first tried) produces zero
  # actual refit failures at this seed -- rchisq() draws from a near-zero
  # ncp are still almost surely positive, so .ppeChiSquareMle() almost
  # always has enough retained values to fit. That made
  # `n_successful + n_failed == bootstrapSamples` tautological: it holds by
  # construction regardless of whether the tryCatch() error handling in
  # .ppeParametricBootstrap() does anything at all. df = 1000 with only 2
  # retained values instead reliably starves most refits (fewer retained
  # values than the density has effective degrees of freedom), producing
  # real, non-zero failures to count.
  res <- .ppeParametricBootstrap(1e10, df = 1000, nRetained = 2L,
                                 bootstrapSamples = 30L, seed = 5L)

  expect_equal(res$n_successful + res$n_failed, 30L)
  expect_gt(res$n_failed, 0L)
  expect_gt(res$n_successful, 0L)
})

test_that("the df-target bootstrap also brackets its point estimate", {
  # Only the ncp-target path is exercised above (that is what the plan's
  # tests cover); a Type-I comparison bootstraps df instead (ncp held at 0),
  # a separate code path in .ppeParametricBootstrap() that needs its own
  # coverage.
  res <- .ppeParametricBootstrap(4, df = NULL, nRetained = 300L,
                                 bootstrapSamples = 200L, target = "df",
                                 seed = 9L)

  expect_lt(res$ci_lower, 4)
  expect_gt(res$ci_upper, 4)
  expect_equal(res$interval_type, "model_based")
})

test_that(".ppeDefaultSeed wraps correctly for a seed near the integer maximum", {
  # base + offset done in integer arithmetic would silently overflow to NA
  # (with a warning) once runInfo$seed sits close to
  # .Machine$integer.max (2147483647L); .ppeDefaultSeed() must do the
  # addition in double precision and only cast back to integer once reduced
  # modulo 2147483647.
  sse <- fake_ppe_sse_object(df = 1, ncp = 8, n = 20L, seed = 2147483000L)
  cmp <- sseComparison("simulation", "alt1", df = 1)

  expect_silent(seed <- .ppeDefaultSeed(sse, cmp))
  expect_false(is.na(seed))
  expect_true(is.integer(seed))
})

test_that(".ppeDefaultSeed distinguishes permutation-related labels", {
  # sum(utf8ToInt(label)) is permutation-invariant, so two labels that are
  # anagrams of each other -- not a contrived case, but a natural pair for
  # the two directions of one comparison -- collided on the same offset
  # (887 for both "dose A vs B" and "dose B vs A"). Position-weighting each
  # character's code by its 1-based index breaks that symmetry.
  sse <- fake_ppe_sse_object(df = 1, ncp = 8, n = 20L, seed = 42L)
  cmp_ab <- sseComparison("simulation", "alt1", df = 1, label = "dose A vs B")
  cmp_ba <- sseComparison("simulation", "alt1", df = 1, label = "dose B vs A")

  expect_false(identical(.ppeDefaultSeed(sse, cmp_ab), .ppeDefaultSeed(sse, cmp_ba)))
})

# --- ppeSummary() ----------------------------------------------------------

test_that("ppeSummary reports estimate, interval, counts, and provenance", {
  sse <- fake_ppe_sse_object(df = 1, ncp = 10, n = 80L, seed = 21L)
  cmp <- sseComparison("simulation", "alt1", df = 1)

  s <- ppeSummary(sse, comparisons = cmp, bootstrapSamples = 100L, bootSeed = 4L)

  expect_equal(s$parameter, "ncp")
  expect_equal(s$mode, "power")
  expect_equal(s$df_source, "explicit")
  expect_equal(s$interval_type, "model_based")
  expect_true(all(c("n", "n_nonpositive", "n_bootstrap_successful",
                    "n_bootstrap_failed", "boundary", "probability") %in% names(s)))
  expect_gt(s$estimate, 0)
})

test_that("a comparison with a custom critical value is refused by PPE", {
  sse <- fake_ppe_sse_object(df = 1, ncp = 10, n = 40L, seed = 22L)
  cmp <- sseComparison("simulation", "alt1", criticalValue = 5)

  expect_error(ppeSummary(sse, comparisons = cmp), "noncentral chi-square")
})

test_that("ppeSummary(bootstrapSamples = 0) keeps the point estimate and counts, drops only the interval", {
  # A user must be able to get the cheap point estimate without paying for
  # the bootstrap.
  sse <- fake_ppe_sse_object(df = 1, ncp = 10, n = 60L, seed = 23L)
  cmp <- sseComparison("simulation", "alt1", df = 1)

  s <- ppeSummary(sse, comparisons = cmp, bootstrapSamples = 0L)

  expect_true(is.na(s$ci_lower))
  expect_true(is.na(s$ci_upper))
  expect_equal(s$n_bootstrap_successful, 0L)
  expect_equal(s$n_bootstrap_failed, 0L)
  expect_gt(s$estimate, 0)
})

test_that("ppeSummary reports one row per comparison, distinguishing power from type1", {
  power_sse <- fake_ppe_sse_object(df = 1, ncp = 10, n = 60L, seed = 24L)
  type1_sse <- fake_ppe_type1_sse_object(df = 1, n = 60L, seed = 24L)

  power_cmp <- sseComparison("simulation", "alt1", df = 1, label = "power check")
  type1_cmp <- sseComparison("alt1", "simulation", df = 1, label = "type1 check")

  power_row <- ppeSummary(power_sse, comparisons = power_cmp, bootstrapSamples = 30L, bootSeed = 1L)
  type1_row <- ppeSummary(type1_sse, comparisons = type1_cmp, bootstrapSamples = 30L, bootSeed = 1L)

  s <- rbind(power_row, type1_row)

  expect_equal(nrow(s), 2L)
  expect_equal(s$mode, c("power", "type1"))
  expect_equal(s$parameter, c("ncp", "df"))
  expect_true(all(s$interval_type == "model_based"))
})
