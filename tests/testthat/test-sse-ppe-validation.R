skip_on_cran()

test_that("proportional fixtures give a consistent lambda per subject", {
  runs <- list(
    fake_ppe_sse_object(df = 1, ncp = 5,  n = 400L, subjects = 20L, seed = 51L),
    fake_ppe_sse_object(df = 1, ncp = 10, n = 400L, subjects = 40L, seed = 52L),
    fake_ppe_sse_object(df = 1, ncp = 20, n = 400L, subjects = 80L, seed = 53L)
  )
  v <- validateSSEPpeScaling(
    runs, comparisons = sseComparison("simulation", "alt1", df = 1)
  )

  expect_equal(v$table$lambda_per_subject, rep(0.25, 3), tolerance = 0.15)
  expect_false(v$nonlinear)
  # A near-perfect through-origin fit: residual variance around the fitted
  # line is small relative to the total variance in lambda across runs.
  # Numerically verified at ~0.0002 for this fixture -- 0.05 leaves a wide
  # margin while still catching a wrong lack-of-fit formula (e.g. a swapped
  # numerator/denominator, which would push this well past 1).
  expect_lt(v$lackOfFit, 0.05)
})

test_that("a deliberately nonlinear fixture is flagged", {
  runs <- list(
    fake_ppe_sse_object(df = 1, ncp = 5, n = 400L, subjects = 20L, seed = 54L),
    fake_ppe_sse_object(df = 1, ncp = 8, n = 400L, subjects = 40L, seed = 55L),
    fake_ppe_sse_object(df = 1, ncp = 9, n = 400L, subjects = 80L, seed = 56L)
  )
  expect_warning(
    v <- validateSSEPpeScaling(
      runs, comparisons = sseComparison("simulation", "alt1", df = 1)
    ),
    "not proportional"
  )
  expect_true(v$nonlinear)
  # A poor through-origin fit: numerically verified at ~1.78 for this
  # fixture. 0.5 stays well below that while still well above the
  # proportional fixture's ~0.0002, giving real separation between the two.
  expect_gt(v$lackOfFit, 0.5)
})

test_that("two study sizes give only a descriptive ratio", {
  runs <- list(
    fake_ppe_sse_object(df = 1, ncp = 5,  n = 200L, subjects = 20L, seed = 57L),
    fake_ppe_sse_object(df = 1, ncp = 10, n = 200L, subjects = 40L, seed = 58L)
  )
  v <- validateSSEPpeScaling(
    runs, comparisons = sseComparison("simulation", "alt1", df = 1)
  )

  expect_true(is.na(v$lackOfFit))
  expect_equal(nrow(v$table), 2L)
})

test_that("incompatible runs are refused rather than silently pooled", {
  # An explicit `comparisons` argument is shared across every run by design
  # (see the roxygen docs), so it can never disagree with itself -- the only
  # way two runs' resolved comparisons can genuinely differ is when each run
  # supplies its OWN declared comparison via `runInfo$comparisons` and
  # `comparisons = NULL` lets each fall back to its own. That is the real
  # scenario this check guards: someone hands in a batch of runs under the
  # (false) assumption they all mean the same comparison.
  run_a <- fake_ppe_sse_object(df = 1, ncp = 5, n = 200L, subjects = 20L, seed = 59L)
  run_b <- fake_ppe_sse_object(df = 1, ncp = 5, n = 200L, subjects = 40L, seed = 60L)
  run_b$runInfo$comparisons <- sseComparison("simulation", "alt1", df = 2)

  # run_a has no declared comparisons of its own, so it also triggers the
  # unrelated "df inferred from parameter counts" warning on the way to the
  # incompatibility error -- expected and not what this test is checking.
  expect_error(
    suppressWarnings(validateSSEPpeScaling(list(run_a, run_b))),
    "incompatible"
  )
})

# @noRd -- the plan sketch's `if (!is.list(runs) || length(runs) < 2L)` guard
# is not exercised by any of the four tests above, but it mirrors the
# argument validation every other exported PPE function in this package
# performs up front (e.g. .assertSSEObject() gating ppeSummary()); a `runs`
# that isn't a usable list must fail with a clear message, not a cryptic
# internal error several frames down.
test_that("fewer than two runs is refused", {
  expect_error(
    validateSSEPpeScaling(list(fake_ppe_sse_object(df = 1, ncp = 5)),
                          comparisons = sseComparison("simulation", "alt1", df = 1)),
    "at least 2"
  )
  expect_error(
    validateSSEPpeScaling("not a list",
                          comparisons = sseComparison("simulation", "alt1", df = 1)),
    "at least 2"
  )
})

test_that("Type-I comparisons are rejected: proportional scaling is a power-curve concept", {
  # .ppeFit()'s target split: a power comparison estimates ncp (df fixed) --
  # the value this validator scales with study size. A Type-I comparison
  # instead fixes ncp = 0 and estimates df, so there is no noncentrality to
  # scale in the first place; "lambda per subject" would be meaningless for
  # it. Reject explicitly rather than silently reporting a nonsense ratio.
  runs <- list(
    fake_ppe_type1_sse_object(df = 1, n = 200L, seed = 61L, subjects = 20L),
    fake_ppe_type1_sse_object(df = 1, n = 200L, seed = 62L, subjects = 40L)
  )
  expect_error(
    validateSSEPpeScaling(
      runs, comparisons = sseComparison("alt1", "simulation", df = 1)
    ),
    "power"
  )
})

test_that("a run resolving to more than one comparison is refused, not silently truncated", {
  # comparisons = NULL falls back to legacy inference, which returns one
  # comparison per alternative model. fake_ppe_mixed_sse_object() has two
  # (alt1, alt2), so leaving comparisons unset here is ambiguous: silently
  # taking the first would validate a comparison the caller never named.
  runs <- list(
    fake_ppe_mixed_sse_object(df = 1, powerNcp = 5, n = 200L, seed = 63L, subjects = 20L),
    fake_ppe_mixed_sse_object(df = 1, powerNcp = 10, n = 200L, seed = 64L, subjects = 40L)
  )
  expect_error(
    suppressWarnings(validateSSEPpeScaling(runs)),
    "exactly one comparison"
  )
})
