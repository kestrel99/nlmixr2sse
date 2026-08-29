# Distribution-Based PPE Design

## Goal

Add the distribution-based parametric power estimation (PPE) method of Ueckert,
Karlsson and Hooker (2016) — the estimator PsN implements — to `nlmixr2sse` as
the new default, retaining the existing exceedance estimator under an explicit
method name, together with explicit model comparisons, explicit degrees of
freedom, a parametric-bootstrap uncertainty interval, a Type-I/df-estimation
mode, and an ECDF-versus-fitted distribution diagnostic.

This design is executed as Tasks 2 and 5–8 of the integrated plan at
`docs/superpowers/plans/2026-08-30-sse-statistical-rigour.md`, which merges it
with the broader statistical-rigour work. Where the two disagreed, the
reconciliation is recorded under **API** below.

## Background

### What nlmixr2sse does today

`.ppePowerPlotData()` (`R/plot-sse.R:350`) estimates the noncentrality
parameter by **matching a single empirical exceedance probability**:

```r
successes  <- sum(test_stat > threshold)
point_prob <- .clipProbability(successes / total, total)
point_ncp  <- .solveNcpForProbability(point_prob, threshold = threshold, df = df)
interval   <- .binomialInterval(successes, total, conf.level = conf.level)
```

`.solveNcpForProbability()` inverts `pchisq(threshold, df, ncp, lower.tail = FALSE)`
by `uniroot()`. Uncertainty comes from a Clopper-Pearson interval on the
proportion.

This uses one bit of information per threshold — the proportion exceeding it —
and discards the shape of the ΔOFV distribution entirely.

### What PsN does

`PsN/R-scripts/sse_default.R`:

1. **Retains positive ΔOFVs.** `chi_square_est()` filters `dofvs[dofvs > 0]`.
2. **Estimates by maximum likelihood** from the full retained distribution:
   `optim(par = init, fn = function(ncp) -sum(dchisq(dofvs, df, ncp, log = TRUE)),
   lower = 1e-16, method = "L-BFGS-B")`, with `init = mean(dofvs) - df`.
3. **Parametric bootstrap for uncertainty.** `param_boot_ppe()` draws
   `rchisq(nmc_samples, df, ncp)` `n.boot = 1000` times, re-estimates, and takes
   the 2.5%/97.5% quantiles. `nmc_samples` is the retained (positive) count.
4. **ECDF diagnostic.** Plots `stat_ecdf()` of the ΔOFVs against the fitted
   `pchisq(grid, df, ncp)` with a ribbon from the bootstrap CI.

PsN also runs a **Type-I mode**: when the simulation model is the reduced one,
it estimates `df` by MLE holding `ncp = 0` and reports the Type-I error rate
instead of power. It selects this mode from the model *filename* suffix
(`_base|_r|_red|_reduced`).

PsN's significance cutoff is fixed at `qchisq(1 - alpha, df)`; `nlmixr2sse`
facets over user-supplied ΔOFV thresholds.

### Verified properties we inherit

Each of these was reproduced directly before being written down.

**Truncation without renormalisation.** `chi_square_est()` drops ΔOFV ≤ 0 and
then fits the *unconditional* noncentral chi-square density. For a genuine
noncentral chi-square `P(X > 0) = 1`, so the discarded observations are ones the
model says cannot occur; dropping them and fitting unconditionally biases the
ncp estimate upward. PsN's mitigation is to report `#(dOFV<0)` in its table.

**Boundary solutions are common at low ncp.** When the retained sample mean
falls below `df`, the method-of-moments start is negative and the constrained
MLE sits at the `1e-16` boundary, so estimated power equals `alpha` exactly.
Measured over 200 replicates per case:

| true ncp | df | n | replicates with a boundary solution |
| --- | --- | --- | --- |
| 1 | 5 | 50 | 0 / 200 |
| 2 | 8 | 40 | 1 / 200 |
| 3 | 10 | 30 | 1 / 200 |
| 0.5 | 4 | 60 | **31 / 200** |

This is the constrained MLE behaving correctly, not a defect — but a flat curve
at `power = alpha` and a degenerate interval must be visibly flagged rather than
silently plotted.

**The out-of-bounds start is harmless.** `init = mean(dofvs) - df` can be
negative, below `lower = 1e-16`. `optim()` silently projects it into the
feasible region and converges normally (verified: `init = -3` on `(p-5)^2`
returns `par = 5`, `convergence = 0`). No guard is needed, and adding one does
not change any estimate.

**The bootstrap does not reproduce the truncation.** `param_boot_ppe()` draws
from `rchisq()`, which is strictly positive, so no bootstrap replicate ever
discards an observation — while the real data did. The interval therefore
reflects sampling variability of the estimator under the fitted model, not the
data-generating reality, and understates uncertainty when negative ΔOFVs are
common.

### The sign-convention trap

`.ofvDeltaPlotData()` computes `delta_ofv = OFV(simulation) - OFV(alternative)`
(`R/plot-sse.R:203`), and `.ppePowerPlotData()` negates it, so the test
statistic is `OFV(alternative) - OFV(simulation)`.

That is ≥ 0 in expectation **only when the simulation model is the full model**.
In a Type-I run the simulation model is the reduced one, the sign flips, every
ΔOFV becomes negative, truncation discards all of them, and the estimator has
nothing to fit. Explicit pairing is what makes Type-I mode possible at all.

### The df-inference fragility

`.modelDegreesFreedom()` (`R/plot-sse.R:265`) infers df from schema parameter
counts and has two silent failure paths: it returns `1L` when the model label is
not found in `fitSpecs`, and it clamps with `max(df, 1L)`, so a comparison whose
alternative has *more* parameters than the simulation model also becomes `df = 1`.
Neither emits any signal.

## Chosen Approach

Extract a PPE estimation core into a new `R/ppe.R` containing pure statistics
with no `ggplot2` dependency. `R/plot-sse.R` (currently 815 lines) consumes it
for rendering. The statistics are the risky part of this feature and become
testable without constructing a plot.

The legacy `exceedance` estimator is retained in full under an explicit method
name, for backward compatibility and sensitivity analysis, but is no longer the
default. Its outputs gain unambiguous names
(`threshold_exceedance_probability`, `threshold_implied_ncp`) so they cannot be
confused with the distribution fit.

## API

```r
plotSSEPpePower(
  x,
  comparisons = NULL,
  method      = c("distribution_mle", "exceedance"),
  df          = NULL,
  alpha       = 0.05,
  thresholds  = NULL,
  models      = NULL,
  studySizes  = NULL,
  targetPower = 99,
  conf.level  = 0.95,
  nonpositive = c("warn", "error", "drop"),
  bootstrapSamples = 1000L,
  bootSeed    = NULL,
  diagnostics = FALSE,
  ...
)

plotSSEPpeDiagnostics(x, comparisons = NULL, df = NULL, alpha = 0.05,
                      conf.level = 0.95, bootstrapSamples = 1000L,
                      bootSeed = NULL, ...)

ppeSummary(x, comparisons = NULL, method = c("distribution_mle", "exceedance"),
           df = NULL, alpha = 0.05, conf.level = 0.95,
           nonpositive = c("warn", "error", "drop"),
           bootstrapSamples = 1000L, bootSeed = NULL, ...)

sseComparison(full, reduced, df = NULL, alpha = 0.05,
              criticalValue = NULL, label = NULL)
```

`method` lists `distribution_mle` first so `match.arg()` makes it the default.
This is a **breaking change to plot output**: existing `plotSSEPpePower()` calls
switch estimator and their numbers change. It requires a NEWS entry.

This reconciles two earlier decisions that disagreed. The
2026-08-29 statistical-rigour plan specified `distribution_mle` as the
documented default; an intermediate decision during this design kept the legacy
estimator as the default. The rigour plan wins: the MLE uses the whole ΔOFV
distribution rather than a single exceedance count, and shipping the weaker
estimator by default would be indefensible once the better one exists.

Method names follow the rigour plan (`distribution_mle`, `exceedance`) rather
than PsN-derived names. A name like `psn_mle` would overclaim, since this
package does not seek PsN execution, seed, or file-format parity — only the
same estimator, which is Ueckert, Karlsson and Hooker (2016).

`bootstrapSamples`, `bootSeed`, `alpha`, and `nonpositive` apply only to
`distribution_mle`. Supplying them with `exceedance` is an error rather than a
silent no-op. `bootstrapSamples = 0` skips the bootstrap without changing the
point estimate. `conf.level`
applies to both, selecting the Clopper-Pearson level under
`exceedance` and the bootstrap percentile quantiles under `distribution_mle`
(the `0.95` default gives PsN's 2.5%/97.5%).

`plotSSEPpeDiagnostics()` has no `method` argument. The diagnostic asks whether
the fitted noncentral chi-square describes the observed ΔOFV distribution,
which is only a meaningful question for the MLE fit; it always uses `distribution_mle`.

`df` accepts a scalar, applied to every comparison, or a vector named by
comparison label, matched by name. An unnamed vector of length greater than one
is an error, since positional matching against comparisons would be silent and
order-dependent.

`models` and `comparisons` are mutually exclusive. `models` filters the default
pairing when `comparisons = NULL`; supplying both is an error, because
`comparisons` already states the full set of pairs explicitly.

## Component Design

| Function | File | Responsibility |
| --- | --- | --- |
| `sseComparison()` | `R/ppe.R` | Construct one explicit `full`/`reduced` pair with optional `df` and `label`. |
| `.resolveComparisons()` | `R/ppe.R` | Return the comparison list; when `NULL`, reproduce today's default (each alternative vs the simulation model as `full`). |
| `.comparisonMode()` | `R/ppe.R` | Derive `"power"` or `"type1"` from which member is the simulation model; abort when neither is. |
| `.resolveComparisonDf()` | `R/ppe.R` | Apply the df precedence chain and record `dfSource`; warn on each silent inference path. |
| `.comparisonDeltaOfv()` | `R/ppe.R` | Build `OFV(reduced) - OFV(full)` per sample for one comparison. |
| `.ppeChiSquareMle()` | `R/ppe.R` | PsN's estimator: truncate at 0, maximise the unconditional log-likelihood, estimate `ncp` or `df`. |
| `.ppeParametricBootstrap()` | `R/ppe.R` | Percentile CI by re-estimating from `rchisq()` draws. |
| `.ppeFit()` | `R/ppe.R` | Orchestrate one comparison into a `ppeFit` record. |
| `ppeSummary()` | `R/ppe.R` | One row per comparison, the public reporting surface. |
| `.ppePowerPlotData()` | `R/plot-sse.R` | Gain a `method` branch; `exceedance` path unchanged. |
| `plotSSEPpePower()` | `R/plot-sse.R` | Render power curves, or Type-I point-ranges, optionally returning both plots. |
| `plotSSEPpeDiagnostics()` | `R/plot-sse.R` | ECDF versus fitted CDF with bootstrap ribbon. |

### The `ppeFit` record

```
comparison   label
mode         "power" | "type1"
df           numeric
dfSource     "explicit" | "argument" | "inferred"
parameter    "ncp" | "df"
estimate     numeric
ciLower      numeric
ciUpper      numeric
n            integer   total ΔOFVs
nRetained    integer   ΔOFV > 0
nNegative    integer   n - nRetained
alpha        numeric
threshold    numeric
boundary     logical   estimate pinned at the lower bound
fittedCdf    function(q) -> cumulative probability
```

### Mode derivation

| Simulation model is | Mode | Estimated | Reported |
| --- | --- | --- | --- |
| the `full` member | `power` | `ncp` (df held) | power at `threshold` |
| the `reduced` member | `type1` | `df` (ncp held at 0) | Type-I rate at `alpha` |
| neither member | — | abort | — |

Deriving the mode structurally replaces PsN's filename-suffix heuristic, which
does not port to labelled `nlmixr2sse` models.

### df precedence

1. `df` supplied to `sseComparison()` — `dfSource = "explicit"`.
2. The `df` argument, scalar or named by comparison label — `dfSource = "argument"`.
3. Schema-count inference via `.modelDegreesFreedom()` — `dfSource = "inferred"`.

Inference warns when the label is missing from `fitSpecs` (fallback to `1L`) or
when `max(df, 1L)` clamps a non-positive difference, naming the comparison and
the value used.

### Threshold resolution

| method | `thresholds` supplied | default |
| --- | --- | --- |
| `exceedance` | honoured | positive thresholds from `powerSummary` (unchanged) |
| `distribution_mle` | honoured | `qchisq(1 - alpha, df)`, per comparison |

Because ncp is estimated from the ΔOFV distribution rather than from an
exceedance count, it does not depend on the threshold. Faceting over several
thresholds therefore remains meaningful under `distribution_mle`: one fit, several
readouts.

### Rendering

**Power comparisons.** Unchanged scaling `ncp_N = ncp * N / N_base`, evaluated as
`pchisq(threshold, df, ncp_N, lower.tail = FALSE)`. Ribbon from the bootstrap CI
rather than the binomial interval.

**Type-I comparisons.** No sample-size curve exists (PsN sets
`max_subjects = NA`). Rendered as a point-range of the estimated Type-I rate with
its bootstrap CI per comparison, against a dashed nominal-`alpha` reference line.

A mixed set of power and Type-I comparisons renders the power panel and warns
that Type-I comparisons were omitted from it, directing the reader to
`ppeSummary()`.

**Diagnostic.** `stat_ecdf()` of the retained ΔOFVs over the fitted CDF, with a
ribbon spanning the CI-implied CDFs, faceted by comparison. The subtitle carries
`n`, `nNegative`, `df`, and `dfSource`.

### Reproducibility

The bootstrap consumes random numbers, so an unseeded plot would draw a
different ribbon on every call. `bootSeed` defaults to a value derived from
`x$runInfo$seed`, making a given run's plot deterministic while remaining
distinct across runs. An explicit `bootSeed` overrides it. The RNG state is
restored on exit so plotting never perturbs a caller's stream.

### Honesty surface

The three inherited properties are reported rather than buried:

- `nNegative` appears as a `ppeSummary()` column and in the diagnostic subtitle.
- `boundary = TRUE` emits a warning naming the comparison, and annotates the
  panel, because power then equals `alpha` and the interval is degenerate.
- The bootstrap's failure to reproduce truncation is documented under a
  "Limitations" heading in `?plotSSEPpePower`, alongside the upward bias from
  fitting the unconditional density to truncated data.

## Error Handling

| Condition | Behaviour |
| --- | --- |
| Neither comparison member is the simulation model | Abort naming the comparison; no hypothesis is known true. |
| A comparison names an unknown model label | Abort listing the available labels. |
| `full` and `reduced` are the same model | Abort. |
| `bootstrapSamples`, `bootSeed`, or `alpha` given with `exceedance` | Abort; they have no meaning for that estimator. |
| Both `models` and `comparisons` supplied | Abort; `comparisons` already names every pair. |
| `df` given as an unnamed vector of length > 1 | Abort; positional matching would be silent and order-dependent. |
| Fewer than 2 retained (positive) ΔOFVs | Abort naming the comparison and the retained count; the MLE is not identifiable. |
| All ΔOFVs non-positive | Abort, and state the sign-convention trap explicitly as the likely cause. |
| `optim()` fails to converge | Abort naming the comparison and the optimiser message. |
| Estimate at the lower bound | Warn, set `boundary = TRUE`, continue. |
| df inferred via a silent fallback or clamp | Warn naming the comparison and the df used. |
| `df` argument named for a label with no comparison | Abort naming the unmatched labels. |

## Testing

Unit tests for the estimation core, no `ggplot2` required:

- `.ppeChiSquareMle()` recovers a known `ncp` from `rchisq(n, df, ncp)` within
  tolerance, across several `df`.
- It recovers a known `df` with `ncp` held at 0.
- It reproduces PsN's numbers exactly on a fixed ΔOFV vector, checked against
  values computed from the `sse_default.R` algorithm transcribed into the test.
- Truncation: negative ΔOFVs are excluded, and `nNegative` counts them.
- A sample whose retained mean falls below `df` yields `boundary = TRUE` and a
  warning.
- `.ppeParametricBootstrap()` is deterministic under a fixed seed, and its
  interval brackets the point estimate.
- The RNG stream is restored after a bootstrap.

Comparison and df layer:

- `comparisons = NULL` reproduces today's pairing exactly.
- Mode derivation returns `power`, `type1`, and aborts for the third case.
- ΔOFV is `OFV(reduced) - OFV(full)` under both modes, so a Type-I pairing keeps
  ΔOFV ≥ 0 in expectation — asserted directly against the sign-flip failure.
- df precedence: explicit beats argument beats inferred.
- Each inference warning fires, naming the comparison.

Rendering:

- `plotSSEPpePower(method = "distribution_mle")` returns a `ggplot`.
- `diagnostics = TRUE` returns a two-element named list of `ggplot`s.
- A Type-I comparison renders a point-range, not a curve.
- `plotSSEPpeDiagnostics()` carries an ECDF layer.
- `ppeSummary()` returns one row per comparison with the documented columns.

Regression guard:

- Every existing `test-plot-sse.R` assertion passes untouched, and
  `.ppePowerPlotData(method = "exceedance")` is byte-identical to the
  current output for the same inputs.

Deliberately **not** tested: that `distribution_mle` and `exceedance` agree. They
are different estimators using different information and will not agree except
asymptotically; asserting otherwise would encode a false expectation.

## Risks

- **Inherited statistical bias.** The truncation and bootstrap limitations are
  PsN's, retained deliberately for parity. The mitigation is disclosure, not
  correction; silently "improving" the estimator would defeat the purpose of a
  PsN-compatible mode.
- **Boundary solutions at low power.** Common enough (15.5% in the measured
  case) that users will meet them. Mitigated by the warning and annotation.
- **Type-I fixtures.** Testing Type-I mode needs an SSE run whose simulation
  model is the reduced one. No such fixture exists; one must be built, and it is
  the largest single piece of test work here.
- **`plot-sse.R` growth.** Even with the core extracted, the rendering additions
  push it toward 1000 lines. If it crosses that during implementation, split the
  PPE rendering into `R/plot-ppe.R` rather than letting it sprawl.

## Out of Scope

- Changing the `exceedance` estimator's algorithm. It is retained verbatim, and
  only loses default status and gains unambiguous output field names.
- Non-linear ncp scaling with study size; `ncp ∝ N` is retained from both
  implementations.
- PsN's PDF report layout and `plot_table()` output.
- Automatic model pairing by filename convention; pairing is explicit or
  defaults to today's simulation-versus-alternative behaviour.
- Correcting the truncation bias by fitting a properly renormalised conditional
  likelihood. That is a genuine improvement over PsN and deserves its own spec
  and validation, not a silent divergence inside a parity mode.
