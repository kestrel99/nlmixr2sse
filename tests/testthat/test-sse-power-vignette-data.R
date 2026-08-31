skip_on_cran()

# testthat::test_path("..", "..", ...) below reaches into the source tree
# (vignettes/, inst/scripts/, README.md) from tests/testthat/'s location.
# That layout only exists when tests run against the source tree itself
# (devtools::test()/load_all()); under R CMD check, tests run against the
# INSTALLED package copy, where only tests/ and inst/'s own contents are
# present at that relative location -- vignettes/ and README.md are never
# installed at all. These are repository-hygiene checks, not package
# behavior, so skip them outright when the source tree isn't there to check.
.sse_power_source_tree_available <- dir.exists(testthat::test_path("..", "..", "R"))

test_that("sse power vignette example artifact is available and well-formed", {
  path <- system.file("extdata", "sse_power_example.rds", package = "nlmixr2sse")

  expect_match(path, "sse_power_example[.]rds$")

  obj <- readRDS(path)

  expect_true(is.list(obj))
  expect_equal(
    sort(names(obj)),
    sort(c(
      "effect_plot_data",
      "n_for_effect_plot",
      "power_curve",
      "study_sizes",
      "true_effect"
    ))
  )
  expect_s3_class(obj$power_curve, "data.frame")
  expect_s3_class(obj$effect_plot_data, "data.frame")
  expect_equal(obj$study_sizes, sort(unique(obj$power_curve$study_size)))
  expect_true(all(c("study_size", "power", "n_converged") %in% names(obj$power_curve)))
  expect_true(all(c("study_size", "estimate", "significant") %in% names(obj$effect_plot_data)))
})

test_that("sse power vignette source exists", {
  skip_if_not(
    .sse_power_source_tree_available,
    "source-tree-only check: not meaningful against an installed package copy"
  )
  path <- testthat::test_path("..", "..", "vignettes", "sse-power.Rmd")
  expect_true(file.exists(path))
})

test_that("README mentions the new power vignette", {
  skip_if_not(
    .sse_power_source_tree_available,
    "source-tree-only check: not meaningful against an installed package copy"
  )
  readme <- readLines(
    testthat::test_path("..", "..", "README.md"),
    warn = FALSE
  )

  expect_true(any(grepl('vignette\\("runSSE"', readme)))
  expect_true(any(grepl('vignette\\("sse-power"', readme)))
})

test_that("standalone SSE power example script exists", {
  skip_if_not(
    .sse_power_source_tree_available,
    "source-tree-only check: not meaningful against an installed package copy"
  )
  path <- testthat::test_path(
    "..",
    "..",
    "inst",
    "scripts",
    "generate-sse-power-example.R"
  )

  expect_equal(file.exists(path), TRUE)
})
