test_that("sseComparison requires distinct labels and one reference definition", {
  expect_error(sseComparison("a", "a", df = 1), "distinct")
  expect_error(sseComparison("a", "b"), "exactly one")
  expect_error(sseComparison("a", "b", df = 1, criticalValue = 3.84), "exactly one")
})

test_that("sseComparison derives the chi-square critical value", {
  cmp <- sseComparison("simulation", "no_cov", df = 1, alpha = 0.05)

  expect_s3_class(cmp, "sseComparison")
  expect_equal(cmp$criticalValue, stats::qchisq(0.95, df = 1))
  expect_equal(cmp$dfSource, "explicit")
  expect_equal(cmp$label, "simulation vs. no_cov")
})

test_that("an explicit criticalValue does not assert PPE validity", {
  cmp <- sseComparison("a", "b", criticalValue = 5)

  expect_equal(cmp$criticalValue, 5)
  expect_null(cmp$df)
  expect_false(cmp$ppeEligible)
})

test_that("the simulation role token resolves to the fitted label", {
  sse <- fake_sse_object()
  cmp <- .resolveComparison(sse, sseComparison("simulation", "alt1", df = 1))

  expect_equal(cmp$full, "fake_sse_fit")
  expect_equal(cmp$mode, "power")
})

test_that("a reduced simulation model gives type1 mode", {
  sse <- fake_sse_object()
  cmp <- .resolveComparison(sse, sseComparison("alt1", "simulation", df = 1))

  expect_equal(cmp$mode, "type1")
})

test_that("an unknown label lists the available labels", {
  sse <- fake_sse_object()
  expect_error(
    .resolveComparison(sse, sseComparison("simulation", "nope", df = 1)),
    "fake_sse_fit"
  )
})

# The plan's original version of this test compared "alt1" against a second
# unknown label ("alt1x") -- but .resolveComparison() checks label existence
# before it interprets which member is the simulation model, so an unknown
# label always reports itself as unknown (see the test above), never as
# "neither member simulated". fake_sse_object() only ever carries one
# simulation model and one alternative, so it cannot produce two *known*
# labels that both fail to match the simulation model. A minimal
# SSE-like list with three fit specs (one simulation, two alternatives)
# is built directly here instead of extending the shared fixture, so this
# edge case stays isolated from every other test that relies on
# fake_sse_object()'s two-model shape.
test_that("a comparison with neither member simulated is rejected", {
  sse <- list(
    runInfo = list(
      fitSpecs = list(
        list(label = "sim", role = "simulation"),
        list(label = "alt1", role = "alternative"),
        list(label = "alt2", role = "alternative")
      )
    )
  )
  expect_error(
    .resolveComparison(sse, sseComparison("alt1", "alt2", df = 1)),
    "neither"
  )
})

test_that("legacy comparisons mark df as inferred and warn once", {
  sse <- fake_sse_object()
  expect_warning(
    cmps <- .resolveComparisons(sse, comparisons = NULL, ppe = TRUE),
    "inferred from parameter counts"
  )
  expect_equal(cmps[[1L]]$dfSource, "parameter_count")
  expect_equal(cmps[[1L]]$full, "fake_sse_fit")
  expect_equal(cmps[[1L]]$reduced, "alt1")
})

test_that("the legacy-comparison warning pluralises the comparison count", {
  sse1 <- fake_sse_object(alt_label = "alt1")
  expect_warning(
    .resolveComparisons(sse1, comparisons = NULL, ppe = TRUE),
    "for 1 comparison\\."
  )

  sse2 <- list(
    runInfo = list(
      fitSpecs = list(
        list(label = "sim", role = "simulation", schema = list()),
        list(label = "alt1", role = "alternative", schema = list()),
        list(label = "alt2", role = "alternative", schema = list())
      )
    )
  )
  expect_warning(
    .resolveComparisons(sse2, comparisons = NULL, ppe = TRUE),
    "for 2 comparisons\\."
  )
})

test_that(".resolveComparisons rejects supplying both comparisons and models", {
  sse <- fake_sse_object()
  expect_error(
    .resolveComparisons(
      sse,
      comparisons = sseComparison("simulation", "alt1", df = 1),
      models = "alt1"
    ),
    "not both"
  )
})

test_that(".resolveComparisons rejects duplicate comparison labels", {
  sse <- fake_sse_object()
  cmp1 <- sseComparison("simulation", "alt1", df = 1, label = "dup")
  cmp2 <- sseComparison("alt1", "simulation", df = 1, label = "dup")
  expect_error(
    .resolveComparisons(sse, comparisons = list(cmp1, cmp2)),
    "unique"
  )
})
