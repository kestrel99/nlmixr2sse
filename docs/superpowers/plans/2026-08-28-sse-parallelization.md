# SSE Parallelization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Expose `rxThreads` through `runSSEControl()` and thread it into the `nlmixr2utils` 0.3 worker-plan API, so parallel SSE runs work at any worker count instead of aborting.

**Architecture:** `runSSE()` resolves `rxThreads` exactly once via `nlmixr2utils::resolveRxThreads()`, then passes that single resolved integer to both `.withWorkerPlan()` call sites, both `.plap()` calls, the run banner, and `run_info`. Because rxode2's thread count changes simulated data even under a fixed seed, the resolved integer is persisted and validated on resume.

**Tech Stack:** R package, `nlmixr2utils` (>= 0.3) worker helpers, `rxode2`, testthat 3e, checkmate, cli.

---

## Background For The Implementer

You do not need prior context on this package. Three facts drive every task:

1. **`nlmixr2utils` 0.3 made a breaking change.** `.withWorkerPlan(workers, expr)` now aborts when `workers x rxode2-threads-per-worker` exceeds the machine's core count. `nlmixr2sse` never passes `rxThreads`, so on a 32-core machine `workers = 4` aborts outright. `runSSEControl()` has no `rxThreads` argument, so users cannot fix it.

2. **Thread count changes results.** rxode2 spreads subjects across threads with per-thread RNG streams. Same model, same seed, different thread count produces different simulated data (verified: max abs difference 0.18 between 1 and 4 threads). So `rxThreads` is part of a run's reproducibility contract, not just a speed knob.

3. **Resolve once, pass the integer.** `"auto"` resolves against the *current* ambient `future` plan. Resolving separately at each call site could yield different values mid-run. Resolve once at the top of `runSSE()` and pass the resulting integer everywhere.

Relevant `nlmixr2utils` exports (all already exported, safe to call with `::`):

- `resolveRxThreads(workers, rxThreads)` — returns an integer.
- `.validateRxThreads(rxThreads)` — accepts `NULL`, `"auto"`, or a positive integer; aborts otherwise.
- `.withWorkerPlan(workers, expr, rxThreads = NULL)`
- `.plap(X, FUN, ..., rxThreads = NULL, .label = NULL)` — note `rxThreads` comes after `...`, so it **must** be passed by name.

## File Structure

- `R/sse-control.R` — `runSSEControl()`. Gains the `rxThreads` argument and its validation.
- `R/run-sse.R` — `runSSE()`. Resolves `rxThreads` once, threads it through both parallel sections, records it in `run_info`, passes it to resume validation.
- `R/sse-helpers.R` — `.workerDescription()`. Gains the resolved thread count for the banner.
- `R/recompute-sse.R` — `.validateResumeRequest()`. Gains the `rxThreads` mismatch check.
- `DESCRIPTION` — `nlmixr2utils` minimum raised to `>= 0.3`.
- `NEWS.md` — records the new option and the reproducibility change.
- `tests/testthat/test-control.R` — control-argument tests.
- `tests/testthat/test-run-sse.R` — resume-validation tests.

## Conventions To Follow

- Tests start with `skip_on_cran()` (already at the top of both test files — do not add a second one).
- Errors are raised with the package-local `.abortSSE()`, which takes a cli-formatted string.
- Existing tests use the helper `capture_sse_error(expr)` from `tests/testthat/helper-fake-fit.R`, which returns the condition object or `NULL`.
- Run a single test file with:
  `Rscript -e "devtools::load_all(quiet=TRUE); testthat::test_file('tests/testthat/<file>', reporter='summary')"`
  Prefix with `NOT_CRAN=true` so `skip_on_cran()` does not skip everything.

---

### Task 1: Add `rxThreads` to `runSSEControl()`

**Files:**
- Modify: `R/sse-control.R`
- Test: `tests/testthat/test-control.R`

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-control.R`:

```r
test_that("runSSEControl exposes rxThreads with an auto default", {
  ctl <- runSSEControl()
  expect_equal(ctl$rxThreads, "auto")

  expect_equal(runSSEControl(rxThreads = 1L)$rxThreads, 1L)
  expect_equal(runSSEControl(rxThreads = 8L)$rxThreads, 8L)
  expect_null(runSSEControl(rxThreads = NULL)$rxThreads)
})

test_that("runSSEControl rejects invalid rxThreads", {
  for (bad in list(0L, -1L, 2.5, "many", c(1L, 2L), NA_integer_)) {
    err <- capture_sse_error(runSSEControl(rxThreads = bad))
    expect_s3_class(err, "error")
    expect_match(conditionMessage(err), "rxThreads")
  }
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `NOT_CRAN=true Rscript -e "devtools::load_all(quiet=TRUE); testthat::test_file('tests/testthat/test-control.R', reporter='summary')"`
Expected: FAIL — `ctl$rxThreads` is `NULL`, so `expect_equal(ctl$rxThreads, "auto")` fails.

- [ ] **Step 3: Add the argument**

In `R/sse-control.R`, add the formal immediately after `workers = NULL,`:

```r
  workers = NULL,
  rxThreads = "auto",
```

Add validation immediately after the existing `nlmixr2utils::.validateWorkers(workers)` line:

```r
  nlmixr2utils::.validateWorkers(workers) # nolint: object_usage_linter.
  nlmixr2utils::.validateRxThreads(rxThreads) # nolint: object_usage_linter.
```

Add the field to the returned `structure()` list, immediately after `workers = workers,`:

```r
      workers = workers,
      rxThreads = rxThreads,
```

Add the roxygen `@param` immediately after the existing `@param workers` line:

```r
#' @param rxThreads Number of rxode2 threads each worker may use. `"auto"`
#'   (the default) divides the machine's cores among the workers. A positive
#'   integer sets the count explicitly; `NULL` defers to rxode2's own default.
#'   Note that rxode2's thread count changes simulated values even under a
#'   fixed seed, so this setting affects reproducibility, not only speed.
```

- [ ] **Step 4: Run test to verify it passes**

Run: `NOT_CRAN=true Rscript -e "devtools::load_all(quiet=TRUE); testthat::test_file('tests/testthat/test-control.R', reporter='summary')"`
Expected: PASS

- [ ] **Step 5: Regenerate docs**

Run: `Rscript -e "devtools::document()"`
Expected: `man/runSSEControl.Rd` updated with the new parameter, exit status 0.

- [ ] **Step 6: Commit**

```bash
git add R/sse-control.R tests/testthat/test-control.R man/runSSEControl.Rd
git commit -m "feat: add rxThreads option to runSSEControl"
```

---

### Task 2: Report the resolved thread count in the run banner

**Files:**
- Modify: `R/sse-helpers.R:99-117`
- Test: `tests/testthat/test-control.R`

`.workerDescription()` builds the `Execution` line of the SSE banner. It should
state the thread count too, so the banner records the full parallel setup.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-control.R`:

```r
test_that("workerDescription reports the resolved rxode2 thread count", {
  expect_match(
    nlmixr2sse:::.workerDescription(1L, rxThreads = 8L),
    "sequential \\(workers = 1\\), 8 rxode2 thread"
  )
  expect_match(
    nlmixr2sse:::.workerDescription(4L, rxThreads = 2L),
    "parallel \\(4 workers\\), 2 rxode2 thread"
  )
  # omitted thread count keeps the bare description
  expect_equal(
    nlmixr2sse:::.workerDescription(1L),
    "sequential (workers = 1)"
  )
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `NOT_CRAN=true Rscript -e "devtools::load_all(quiet=TRUE); testthat::test_file('tests/testthat/test-control.R', reporter='summary')"`
Expected: FAIL — `.workerDescription()` takes only one argument, so the call errors with `unused argument (rxThreads = 8L)`.

- [ ] **Step 3: Add the argument**

In `R/sse-helpers.R`, replace the whole `.workerDescription` function (lines 99-117) with:

```r
.workerDescription <- function(workers, rxThreads = NULL) {
  base <- if (is.null(workers)) {
    "sequential (using current future::plan())"
  } else if (identical(workers, "auto")) {
    if (requireNamespace("future", quietly = TRUE)) {
      paste0(
        "parallel, auto (",
        future::availableCores(omit = 1L),
        " workers)"
      )
    } else {
      "sequential (future not available)"
    }
  } else if (identical(as.integer(workers), 1L)) {
    "sequential (workers = 1)"
  } else {
    paste0("parallel (", workers, " workers)")
  }

  if (is.null(rxThreads)) {
    return(base)
  }
  paste0(
    base,
    ", ",
    rxThreads,
    " rxode2 thread",
    if (identical(as.integer(rxThreads), 1L)) "" else "s",
    "/worker"
  )
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `NOT_CRAN=true Rscript -e "devtools::load_all(quiet=TRUE); testthat::test_file('tests/testthat/test-control.R', reporter='summary')"`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add R/sse-helpers.R tests/testthat/test-control.R
git commit -m "feat: report resolved rxode2 thread count in SSE banner"
```

---

### Task 3: Resolve `rxThreads` once and thread it through `runSSE()`

**Files:**
- Modify: `R/run-sse.R:152-153`, `R/run-sse.R:294` (run_info), `R/run-sse.R:412`, `R/run-sse.R:432`, `R/run-sse.R:480`, `R/run-sse.R:506`

This task has no new test of its own — Task 5's resume tests and the existing
integration tests cover it. Verification is by running the full suite in Step 6.

- [ ] **Step 1: Resolve once, before the banner**

In `R/run-sse.R`, find:

```r
  sim_schema <- .rawResultsSchemaForFit(fit)
  worker_desc <- .workerDescription(control$workers)
```

Replace with:

```r
  sim_schema <- .rawResultsSchemaForFit(fit)
  resolved_rx_threads <- nlmixr2utils::resolveRxThreads(
    control$workers,
    control$rxThreads
  )
  worker_desc <- .workerDescription(control$workers, resolved_rx_threads)
```

- [ ] **Step 2: Record the resolved value in `run_info`**

Find this line in the `run_info` assembly:

```r
  run_info$workers <- control$workers
```

Replace with:

```r
  run_info$workers <- control$workers
  run_info$rxThreads <- resolved_rx_threads
```

Record the **resolved integer**, never the literal `"auto"` — a completed run
must state the thread count it actually used.

- [ ] **Step 3: Pass it to the simulation worker plan**

Find the simulation block opening:

```r
      nlmixr2utils::.withWorkerPlan(control$workers, {
        nlmixr2utils::.plap(
          pending_data,
```

Replace those three lines with:

```r
      nlmixr2utils::.withWorkerPlan(control$workers, {
        nlmixr2utils::.plap(
          pending_data,
          rxThreads = resolved_rx_threads,
```

Then find that block's closing lines:

```r
          .label = function(sample_id) {
            paste0("simulation ", sample_id, "/", samples)
          }
        )
      })
```

Replace with:

```r
          .label = function(sample_id) {
            paste0("simulation ", sample_id, "/", samples)
          }
        )
      }, rxThreads = resolved_rx_threads)
```

- [ ] **Step 4: Pass it to the fitting worker plan**

Find the fitting block opening:

```r
    nlmixr2utils::.withWorkerPlan(control$workers, {
      nlmixr2utils::.plap(
        pending_fit,
```

Replace those three lines with:

```r
    nlmixr2utils::.withWorkerPlan(control$workers, {
      nlmixr2utils::.plap(
        pending_fit,
        rxThreads = resolved_rx_threads,
```

Then find that block's closing lines:

```r
        .label = function(key) {
          parts <- strsplit(key, "::", fixed = TRUE)[[1L]]
          paste0(parts[[1L]], " ", parts[[2L]], "/", samples)
        }
      )
    })
```

Replace with:

```r
        .label = function(key) {
          parts <- strsplit(key, "::", fixed = TRUE)[[1L]]
          paste0(parts[[1L]], " ", parts[[2L]], "/", samples)
        }
      )
    }, rxThreads = resolved_rx_threads)
```

- [ ] **Step 5: Verify the package loads**

Run: `Rscript -e "devtools::load_all(quiet=TRUE); cat('loaded ok\n')"`
Expected: `loaded ok`

- [ ] **Step 6: Verify a parallel run no longer aborts**

This is the defect being fixed. Run:

```bash
Rscript -e '
suppressMessages(library(nlmixr2utils))
cores <- parallel::detectCores()
cat("cores:", cores, " rxode2 default threads:", rxode2::getRxThreads(), "\n")
# what nlmixr2sse used to do: no rxThreads
old <- try(nlmixr2utils::.withWorkerPlan(4, { "ok" }), silent = TRUE)
cat("workers=4, no rxThreads ->", if (inherits(old, "try-error")) "ABORT" else "ok", "\n")
# what it does now: resolved rxThreads
res <- nlmixr2utils::resolveRxThreads(4, "auto")
new <- try(nlmixr2utils::.withWorkerPlan(4, { "ok" }, rxThreads = res), silent = TRUE)
cat("workers=4, rxThreads=", res, " ->", if (inherits(new, "try-error")) "ABORT" else "ok", "\n")
'
```

Expected on a machine where `4 x` rxode2's default exceeds the core count: the
first line reports `ABORT`, the second reports `ok`. On a machine with enough
cores both report `ok`; that is fine and still demonstrates the wiring.

- [ ] **Step 7: Run the full suite for regressions**

Run: `NOT_CRAN=true Rscript -e "devtools::load_all(quiet=TRUE); devtools::test()"`
Expected: no *new* failures. Note that this repository has pre-existing
failures unrelated to this change — see "Pre-Existing Failures" at the end of
this plan. Compare against that list; anything else is a regression you caused.

- [ ] **Step 8: Commit**

```bash
git add R/run-sse.R
git commit -m "fix: thread resolved rxThreads through SSE worker plans"
```

---

### Task 4: Raise the `nlmixr2utils` dependency

**Files:**
- Modify: `DESCRIPTION`

`resolveRxThreads()`, `.validateRxThreads()`, and the `rxThreads` arguments are
all 0.3 features. The current floor of `>= 0.2` would let the package install
against a version where the code cannot work.

- [ ] **Step 1: Raise the minimum**

In `DESCRIPTION`, under `Imports:`, change:

```
    nlmixr2utils (>= 0.2),
```

to:

```
    nlmixr2utils (>= 0.3),
```

- [ ] **Step 2: Verify the installed version satisfies it**

Run: `Rscript -e "cat(as.character(packageVersion('nlmixr2utils')), '\n')"`
Expected: `0.3.1` or higher.

- [ ] **Step 3: Commit**

```bash
git add DESCRIPTION
git commit -m "build: require nlmixr2utils >= 0.3 for rxThreads support"
```

---

### Task 5: Validate `rxThreads` on resume

**Files:**
- Modify: `R/recompute-sse.R:109-115` (signature), and the body before the `existing_labels` block
- Modify: `R/run-sse.R` — the `.validateResumeRequest()` call site
- Test: `tests/testthat/test-run-sse.R`

Because thread count changes simulated data, resuming a partial run under a
different `rxThreads` would mix replicates drawn under two configurations.
Abort on mismatch. Directories written before this change carry no recorded
value; warn rather than abort, so existing work is not stranded.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-run-sse.R`:

```r
test_that("validateResumeRequest aborts on an rxThreads mismatch", {
  run_info <- list(
    fitName = "fake_sse_fit",
    samples = 2L,
    parameterSource = "fixed",
    estimateSimulation = TRUE,
    rxThreads = 8L
  )

  err <- capture_sse_error(
    .validateResumeRequest(
      existingRunInfo = run_info,
      fitName = "fake_sse_fit",
      samples = 2L,
      control = runSSEControl(workers = 1L),
      requestedLabels = character(0),
      rxThreads = 2L
    )
  )

  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "rxThreads")
  expect_match(conditionMessage(err), "8")
  expect_match(conditionMessage(err), "2")
})

test_that("validateResumeRequest warns when rxThreads was never recorded", {
  run_info <- list(
    fitName = "fake_sse_fit",
    samples = 2L,
    parameterSource = "fixed",
    estimateSimulation = TRUE
  )

  expect_warning(
    .validateResumeRequest(
      existingRunInfo = run_info,
      fitName = "fake_sse_fit",
      samples = 2L,
      control = runSSEControl(workers = 1L),
      requestedLabels = character(0),
      rxThreads = 4L
    ),
    "rxThreads"
  )
})

test_that("validateResumeRequest accepts a matching rxThreads", {
  run_info <- list(
    fitName = "fake_sse_fit",
    samples = 2L,
    parameterSource = "fixed",
    estimateSimulation = TRUE,
    rxThreads = 4L
  )

  expect_silent(
    .validateResumeRequest(
      existingRunInfo = run_info,
      fitName = "fake_sse_fit",
      samples = 2L,
      control = runSSEControl(workers = 1L),
      requestedLabels = character(0),
      rxThreads = 4L
    )
  )
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `NOT_CRAN=true Rscript -e "devtools::load_all(quiet=TRUE); testthat::test_file('tests/testthat/test-run-sse.R', reporter='summary')"`
Expected: the three new tests FAIL with `unused argument (rxThreads = ...)`.

- [ ] **Step 3: Add the parameter and the check**

In `R/recompute-sse.R`, change the signature:

```r
.validateResumeRequest <- function(
  existingRunInfo,
  fitName,
  samples,
  control,
  requestedLabels
) {
```

to:

```r
.validateResumeRequest <- function(
  existingRunInfo,
  fitName,
  samples,
  control,
  requestedLabels,
  rxThreads = NULL
) {
```

Then, immediately **before** the line `existing_labels <- .existingModelLabels(existingRunInfo)`, insert:

```r
  if (!is.null(rxThreads)) {
    recorded <- existingRunInfo$rxThreads
    if (is.null(recorded)) {
      cli::cli_warn(c(
        "!" = "This run directory predates {.arg rxThreads} tracking, so the thread count it used is unknown.",
        "i" = "rxode2 thread count changes simulated values, so resumed replicates may not be comparable to the existing ones."
      ))
    } else if (!identical(as.integer(recorded), as.integer(rxThreads))) {
      .abortSSE(
        paste0(
          "Existing run directory used {.arg rxThreads = {as.integer(recorded)}}, ",
          "but this run resolves to {.arg rxThreads = {as.integer(rxThreads)}}. ",
          "rxode2 thread count changes simulated values, so the replicates would not be comparable. ",
          "Use {.code runSSEControl(rxThreads = {as.integer(recorded)})} to resume, or {.code restart = TRUE} for a fresh run."
        )
      )
    }
  }

```

- [ ] **Step 4: Pass it from the call site**

In `R/run-sse.R`, find the `.validateResumeRequest()` call and add the argument.
Change:

```r
    .validateResumeRequest(
      existingRunInfo = existing_run_info,
      fitName = fitName,
      samples = samples,
      control = control,
      requestedLabels = .requestedModelLabels(.fitTaskSpecs(
        fit,
        fitName,
        alternatives,
        control
      ))
    )
```

to:

```r
    .validateResumeRequest(
      existingRunInfo = existing_run_info,
      fitName = fitName,
      samples = samples,
      control = control,
      requestedLabels = .requestedModelLabels(.fitTaskSpecs(
        fit,
        fitName,
        alternatives,
        control
      )),
      rxThreads = resolved_rx_threads
    )
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `NOT_CRAN=true Rscript -e "devtools::load_all(quiet=TRUE); testthat::test_file('tests/testthat/test-run-sse.R', reporter='summary')"`
Expected: the three new tests PASS.

Note: the pre-existing `runSSE resumes cached partial runs without duplicating
rows` test resumes a directory whose `run_info` has no `rxThreads`, so it will
now emit the new warning. If that test fails **because of an unexpected
warning**, wrap its `runSSE(...)` call in `suppressWarnings(...)`. Do not
weaken the warning itself — it is the designed behaviour.

- [ ] **Step 6: Run the full suite**

Run: `NOT_CRAN=true Rscript -e "devtools::load_all(quiet=TRUE); devtools::test()"`
Expected: no new failures beyond the pre-existing list below.

- [ ] **Step 7: Commit**

```bash
git add R/recompute-sse.R R/run-sse.R tests/testthat/test-run-sse.R
git commit -m "feat: validate rxThreads when resuming an SSE run"
```

---

### Task 6: Document the change

**Files:**
- Modify: `NEWS.md`
- Modify: `README.md`

- [ ] **Step 1: Add the NEWS entry**

At the top of `NEWS.md`, immediately under the `# nlmixr2sse 0.1` heading, add:

```markdown
* `runSSEControl()` gains `rxThreads`, which sets how many rxode2 threads each
  parallel worker may use. It defaults to `"auto"`, dividing the machine's
  cores among the workers. Parallel SSE runs previously aborted at worker
  counts where the workers' combined rxode2 threads exceeded the core count,
  and the setting needed to avoid that was not reachable through
  `runSSEControl()`. Requires `nlmixr2utils` >= 0.3.
* **Reproducibility note:** rxode2's thread count changes simulated values even
  under a fixed seed, so `rxThreads` affects results, not only speed. The
  resolved thread count is now recorded in `run_info` and checked when a run is
  resumed. Because the default is `"auto"`, which derives from the machine's
  core count, reproducing a study on different hardware requires passing the
  recorded integer explicitly, e.g. `runSSEControl(rxThreads = 16)`. Runs made
  with earlier versions used rxode2's own default and will not reproduce
  byte-for-byte under the new default.
```

- [ ] **Step 2: Document it in the README**

In `README.md`, immediately before the `## Credit where it's due` heading, add:

```markdown
## Parallel execution

`runSSEControl(workers = )` sets how many replicates run in parallel, and
`rxThreads` sets how many rxode2 threads each worker may use. The default,
`rxThreads = "auto"`, divides the machine's cores among the workers.

```r
sse <- runSSE(
  fit,
  samples = 500,
  seed = 42,
  control = runSSEControl(workers = 4, rxThreads = "auto")
)
```

rxode2's thread count changes simulated values even under a fixed seed, so
`rxThreads` is part of a run's reproducibility contract. The resolved count is
stored in `run_info` and checked on resume. To reproduce a study on different
hardware, pass the recorded integer explicitly rather than relying on `"auto"`.
```

- [ ] **Step 3: Verify the README example is valid R**

Run:

```bash
Rscript -e '
devtools::load_all(quiet = TRUE)
ctl <- runSSEControl(workers = 4, rxThreads = "auto")
stopifnot(identical(ctl$rxThreads, "auto"), identical(ctl$workers, 4))
cat("README example constructs a valid control\n")
'
```

Expected: `README example constructs a valid control`

- [ ] **Step 4: Commit**

```bash
git add NEWS.md README.md
git commit -m "docs: document rxThreads and its reproducibility implications"
```

---

## Pre-Existing Failures

`devtools::test()` on this repository fails before any of this plan's changes.
Do not try to fix these here; use the list to tell your regressions apart from
inherited breakage.

1. `test-run-sse.R` — `simulationRecord preserves row identifiers...`: errors
   with `bad matrix specification`. The fit's empty `sigma` is passed to
   `rxSolve()`, which rejects it.
2. `test-run-sse.R` — `runSSE resumes cached partial runs...`: errors with
   `numbers of columns of arguments do not match`. The fake fit's `$` method
   breaks `fit$iniDf`, producing an empty schema.
3. `test-run-sse.R` — `runSSE addModels reuses saved datasets...`: errors with
   `Could not recover the original estimation data from the reference fit.`
4. `test-run-sse.R` — `runSSE supports raw-results and covariance parameter
   sources end to end`: **skips**, masked by a `cli` `.envir` defect in
   `nlmixr2utils:::.abortRawResults()`.
5. `test-sse-power-vignette-data.R` — two failures asserting the README
   mentions `vignette("runSSE"...)` and `sse-power`. The README does not.

After this plan, expect the same set, except that #2's test may additionally
need `suppressWarnings()` as noted in Task 5 Step 5.

## Self-Review

**Spec coverage:**

- `rxThreads` on `runSSEControl` with `"auto"` default — Task 1.
- Resolve once in `runSSE()` — Task 3 Step 1.
- Pass to both `.withWorkerPlan()` and both `.plap()` calls — Task 3 Steps 3-4.
- Record the resolved integer in `run_info` — Task 3 Step 2.
- `.workerDescription()` reports the thread count — Task 2.
- Resume mismatch aborts; missing record warns — Task 5.
- `DESCRIPTION` bump to `>= 0.3` — Task 4.
- NEWS entry for the option and the reproducibility break — Task 6.
- Error-handling table: invalid `rxThreads` (Task 1), oversubscription passed
  through unwrapped from `nlmixr2utils` (no code needed — verified in Task 3
  Step 6), resume mismatch and missing record (Task 5).
- Spec's "deliberately not tested" note — honoured; no task asserts that two
  thread counts agree.

**Placeholder scan:** No TBD/TODO. Every code step shows complete code; every
run step shows the command and expected output.

**Type consistency:** `resolved_rx_threads` is the same integer-valued name in
Task 3 Steps 1-4 and Task 5 Step 4. `.workerDescription(workers, rxThreads)`
matches its Task 3 Step 1 call. `.validateResumeRequest(..., rxThreads = )`
matches between Task 5 Step 3 and Step 4. `control$rxThreads` is only ever read
raw at the single resolution point; the resolved integer is used everywhere
else.
