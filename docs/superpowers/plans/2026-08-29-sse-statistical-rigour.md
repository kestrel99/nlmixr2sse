# SSE Statistical Rigour Improvement Plan

> **Superseded** by `docs/superpowers/plans/2026-08-30-sse-statistical-rigour.md`,
> which merges this roadmap with the PPE design in
> `docs/superpowers/specs/2026-08-30-sse-psn-ppe-design.md` and expands every
> task into executable TDD steps. Kept for provenance; do not execute from here.

**Goal:** Strengthen the statistical definition, diagnostics, and uncertainty
quantification of `nlmixr2sse` without reproducing NONMEM- or PsN-specific
execution machinery.

**Architecture:** Preserve the existing simulation and estimation engine and
the positive-definite `independent_iw` and `joint` parameter-draw modes. Add an
explicit comparison specification, keep simulation truths separate from
estimation starting values, calculate all model-comparison results on paired
evaluable replicates, add Monte Carlo uncertainty to empirical summaries, and
replace threshold-by-threshold PPE as the recommended method with a
full-distribution noncentral-chi-square likelihood fit plus diagnostics and
sample-size-scaling validation. Retain current behavior behind explicit legacy
options where practical.

**Tech stack:** R package, `nlmixr2est`, `rxode2`, `nlmixr2utils`, `ggplot2`,
`checkmate`, `cli`, testthat 3e, and base R statistical functions.

---

## Scope

This plan addresses statistical properties that can change scientific
conclusions:

- what hypothesis comparison is being made;
- which replicate pairs contribute to each estimate;
- how Monte Carlo error is quantified;
- whether estimation starts alter convergence and selection of usable fits;
- how the PPE noncentrality parameter is estimated;
- whether the assumed noncentral chi-square distribution and linear
  sample-size scaling are empirically credible; and
- whether generated parameter distributions behave as intended.

This plan does **not** seek file-format, random-seed, command-line, or execution
parity with PsN. In particular, it does not add NONMEM covariance-file parsing,
natural-coordinate Gaussian OMEGA draws, `$PRIOR` syntax, MSF handling,
control-stream generation, numbered parameter matching, or PsN output layouts.

## Current statistical gaps

### Implicit hypothesis definitions

The package currently constructs each OFV comparison from the simulation model
and one alternative, infers degrees of freedom from schema parameter counts,
and assigns `power` and `type1` labels from the sign of the OFV difference. It
cannot establish nesting, recognize boundary hypotheses, verify comparable
likelihoods, or know the scientifically intended null model.

### Incomplete Monte Carlo accounting

Model failures and output filters can produce different usable replicate sets.
Every comparison must therefore report the number of attempted, paired,
evaluable, successful, and excluded replicates. A percentage without its
denominator and Monte Carlo interval is insufficient.

Parameter summaries expose `n_effective`, but the field named `rse` is limited
to fixed truths and can become negative when the fixed truth is negative.
Monte Carlo uncertainty should instead be computed from the replicate-level
error statistic and must always be nonnegative.

### Estimation initialization is a numerical intervention

`randomEstimationInits = TRUE` applies each raw-results generating vector to
all compatible fits. This does not change the simulated data distribution, but
starting at the generating truth can change convergence, boundary behavior,
runtime, and which fits enter the summaries. Reference and alternative starts
need to be separately selectable and recorded.

### PPE uses threshold-specific inversion

Current PPE estimates a separate noncentrality parameter for each threshold by
inverting the empirical exceedance probability. It does not fit the complete
OFV-difference distribution. Different thresholds can consequently imply
different effect sizes, and the Clopper-Pearson ribbon assesses only the
binomial exceedance estimate, not the noncentral chi-square assumption.

The distribution-based PPE described by Ueckert, Karlsson, and Hooker (2016)
instead treats the OFV differences as draws from a noncentral chi-square
distribution, estimates a single noncentrality parameter for a specified
degrees of freedom, and scales that parameter with study size. That approach
still requires diagnostics; it is not made valid merely by fitting it.

### Parameter-draw approximations need observable diagnostics

The current covariance modes make intentional approximations:

- `independent_iw` preserves the fitted raw OMEGA mean but discards
  THETA--OMEGA and cross-block dependence and over-disperses nonbinding OMEGA
  diagonals;
- `joint` retains local dependence to first order but can shift raw-scale OMEGA
  means after back-transformation; and
- both use unbounded Gaussian THETA draws.

These modes should remain available, but the realized parameter draws need a
standard adequacy summary rather than relying only on warnings.

## Design decisions

1. **Do not emulate an unconstrained natural-scale OMEGA Normal.** Producing
   indefinite covariance matrices is not a statistical improvement.
2. **Do not rename `independent_iw` as NWPRI.** Its degrees of freedom are
   inferred from standard errors and it is a different distributional model.
3. **Prefer empirical joint draws when justified.** Validated bootstrap or SIR
   parameter vectors supplied through `rawres` are preferable to a local
   covariance approximation when asymmetry or non-Gaussian dependence matters.
4. **Keep initialization separate from the data-generating distribution.** A
   starting-value sensitivity analysis is a numerical operating-characteristic
   analysis, not parameter uncertainty.
5. **Define comparisons explicitly.** Parameter-count differences remain a
   convenience fallback, not an asserted truth.
6. **Use paired evaluable replicates.** Every OFV statistic requires finite,
   accepted results for both models in the same replicate.
7. **Make distribution-based PPE the recommended method.** Preserve the
   current exceedance inversion under an explicit method name for backward
   compatibility and sensitivity analysis.
8. **Never hide negative OFV differences.** A distribution fit may exclude
   values outside noncentral chi-square support, but must report their number
   and fraction and warn or abort according to policy.
9. **Distinguish empirical and model-based uncertainty.** Binomial intervals,
   parameter-summary Monte Carlo intervals, and parametric PPE bootstrap bands
   answer different questions and must have different field names.
10. **Do not automate boundary-test reference distributions.** Users must
   supply an appropriate critical value or use empirical null simulation.
11. **Report marginal operating characteristics by default.** Multiple
    comparisons do not acquire family-wise error control merely by being in the
    same SSE run; any multiplicity procedure must be explicit.
12. **Record every statistical choice in `runInfo` and derived plot data.**

## Proposed public API

### Explicit comparison specification

```r
cmp <- sseComparison(
  full = "simulation",
  reduced = "no_covariate",
  df = 1,
  alpha = 0.05,
  label = "covariate effect"
)

sse <- runSSE(
  fit,
  alternativeModels = alternatives,
  comparisons = list(cmp),
  ...
)
```

`sseComparison()` should accept either:

- `df` and `alpha`, from which the ordinary chi-square critical value is
  calculated; or
- an explicit `criticalValue` for empirical operating characteristics when the
  ordinary chi-square reference is inappropriate.

`"simulation"` is a reserved role token resolved to the actual label of the
fitted simulation model. Other values are explicit persisted model labels.

Distribution-based PPE is enabled only for comparisons with an explicit
positive `df` and an ordinary noncentral-chi-square alternative assumption.
Supplying a custom critical value must not silently assert that PPE remains
valid.

### Role-specific starting values

```r
runSSEControl(
  referenceInitials = "model",
  alternativeInitials = "model"
)
```

Each argument initially supports `"model"` and `"simulation"`. The latter
means that compatible free parameters are initialized from that replicate's
generating vector. Defaults remain `"model"`. The existing
`randomEstimationInits` argument maps to `"simulation"` for both roles with a
lifecycle warning, then is deprecated after a transition period.

This API permits a deliberate sensitivity study such as:

```r
stored_starts <- runSSEControl(
  referenceInitials = "model",
  alternativeInitials = "model"
)

truth_starts <- runSSEControl(
  referenceInitials = "simulation",
  alternativeInitials = "simulation"
)
```

The two runs use the same simulated datasets when their seeds and parameter
sources match, isolating the numerical effect of starting values.

### PPE method and diagnostic API

```r
plotSSEPpePower(
  sse,
  comparisons = "covariate effect",
  method = "distribution_mle",
  bootstrapSamples = 1000,
  conf.level = 0.95
)

plotSSEPpeDiagnostics(sse, comparisons = "covariate effect")
```

The legacy method remains available as `method = "exceedance"`.

### Parameter-draw diagnostics

```r
parameterDrawSummary(sse)
plotSSEParameterDraws(sse)
```

These functions describe the realized generating values. They do not certify
that the approximation is a posterior or exact sampling distribution.

---

## Task 1: Freeze the current statistical contract

**Files:**

- Create: `tests/testthat/test-sse-statistical-contract.R`
- Read: `R/sse-helpers.R`
- Read: `R/plot-sse.R`

- [ ] Add a small synthetic `nlmixr2SSE` fixture with known OFVs, missing OFVs,
  filtered rows, fixed truths, varying truths, and negative truths.
- [ ] Characterize the current sign convention
  `delta = OFV_reference - OFV_alternative` and
  `T = OFV_reduced - OFV_full`.
- [ ] Characterize which rows currently enter parameter and OFV summaries.
- [ ] Characterize current threshold-specific PPE inversion at two thresholds.
- [ ] Characterize `randomEstimationInits = FALSE` and `TRUE` separately for
  reference and alternative fit specifications.
- [ ] Run the complete existing test suite and record baseline failures without
  changing them in this task.

**Acceptance:** Tests describe current behavior precisely and fail only when a
later task intentionally changes that behavior.

## Task 2: Add explicit comparison objects

**Files:**

- Create: `R/sse-comparison.R`
- Create: `tests/testthat/test-sse-comparison.R`
- Modify: `R/run-sse.R`
- Modify: `R/recompute-sse.R`
- Modify: `R/sse-helpers.R`
- Modify: `NAMESPACE`

- [ ] Implement `sseComparison(full, reduced, df = NULL, alpha = 0.05,
  criticalValue = NULL, label = NULL)`.
- [ ] Require distinct, nonempty model labels and a unique comparison label.
- [ ] Require exactly one of an ordinary chi-square definition (`df`, `alpha`)
  or an explicit `criticalValue`.
- [ ] Calculate `criticalValue = qchisq(1 - alpha, df)` only for the ordinary
  chi-square case.
- [ ] Validate comparison labels against the persisted fit specifications.
- [ ] Resolve the reserved `"simulation"` role token to the unique fitted
  simulation-model label and give a clear error when that role was not fitted.
- [ ] Add optional `comparisons` to `runSSE()` and persist normalized objects in
  `runInfo`.
- [ ] Allow `recomputeSSE()` to replace comparison definitions without rerunning
  fits.
- [ ] When comparisons are absent, construct legacy reference-versus-alternative
  comparisons, mark `dfSource = "parameter_count"`, and warn when PPE is
  requested rather than at run time.
- [ ] Do not attempt to prove nesting or likelihood comparability.

**Acceptance:** Every new comparison summary identifies full model, reduced
model, sign convention, degrees of freedom or custom critical value, alpha,
and whether the definition was explicit or inferred.

## Task 3: Compute paired empirical operating characteristics

**Files:**

- Create: `R/sse-comparison-summary.R`
- Create: `tests/testthat/test-sse-comparison-summary.R`
- Modify: `R/sse-helpers.R`
- Modify: `R/sse-methods.R`
- Modify: `R/plot-sse.R`

- [ ] Build one row per comparison and replicate by joining full and reduced
  fits on `sample` after applying the same output-filter policy.
- [ ] Define an evaluable replicate as one with finite accepted OFVs for both
  models.
- [ ] Calculate `test_statistic = OFV_reduced - OFV_full` directly; do not infer
  direction from an ambiguous delta label.
- [ ] Report `n_attempted`, `n_full_evaluable`, `n_reduced_evaluable`,
  `n_paired_evaluable`, `n_excluded`, and `n_exceeding`.
- [ ] Report empirical probability, binomial Monte Carlo standard error, and a
  two-sided exact interval at a configurable confidence level.
- [ ] Retain replicate identifiers for auditing which fits were excluded.
- [ ] Keep legacy `ofvSummary` fields for one release, but build new power plots
  from `comparisonSummary`.
- [ ] Warn when the paired evaluable fraction falls below a configurable
  reporting threshold; do not silently impute failed fits.

**Acceptance:** A comparison with 100 attempted replicates, 90 usable full
fits, 85 usable reduced fits, and 80 paired fits reports denominator 80—not 85,
90, or 100—and its confidence interval is reproducible from the reported count.

## Task 4: Separate estimation starts by model role

**Files:**

- Modify: `R/sse-control.R`
- Modify: `R/sse-helpers.R`
- Modify: `R/run-sse.R`
- Modify: `R/recompute-sse.R`
- Modify: `DESCRIPTION`
- Modify: package-level roxygen file
- Create: lifecycle badge assets through `usethis::use_lifecycle()`
- Modify: `tests/testthat/test-control.R`
- Modify: `tests/testthat/test-run-sse.R`

- [ ] Add `referenceInitials` and `alternativeInitials`, each accepting
  `"model"` or `"simulation"`.
- [ ] Set up the package's currently absent lifecycle infrastructure with
  `usethis::use_lifecycle()` and inspect every generated change.
- [ ] Update `.fitTaskRecord()` according to `spec$role`, rather than using one
  switch for every fit.
- [ ] Keep fixed values unchanged unless `updateFix = TRUE`.
- [ ] Record the resolved starting-value source for each role in `runInfo` and
  validate it on resume and `addModels`.
- [ ] Change the old argument default to `lifecycle::deprecated()`. When it is
  explicitly supplied, map `FALSE` to model starts and `TRUE` to simulation
  starts for both roles, call `lifecycle::deprecate_soft()`, and reject
  contradictory combinations with the new arguments.
- [ ] Verify that changing starting-value policy does not change cached
  parameter vectors or simulated datasets.
- [ ] Add a diagnostic summary of convergence, thrown errors, and output-filter
  exclusions by model and starting-value policy.

**Acceptance:** Reference-only, alternative-only, and both-role simulation
starts are independently testable, and two runs differing only in starts have
identical generating values and simulated datasets.

## Task 5: Implement full-distribution PPE estimation

**Files:**

- Create: `R/sse-ppe.R`
- Create: `tests/testthat/test-sse-ppe.R`
- Modify: `R/plot-sse.R`

- [ ] Implement an internal noncentral-chi-square negative log-likelihood with
  fixed positive degrees of freedom and noncentrality constrained to zero or
  greater.
- [ ] Fit one noncentrality parameter per comparison from finite positive test
  statistics.
- [ ] Report zeros and negative values separately as assumption violations.
- [ ] Add `nonpositive = c("warn", "error", "drop")`, defaulting to `"warn"`;
  `"drop"` suppresses the warning but never suppresses the counts.
- [ ] Use a stable bounded or transformed optimizer and return convergence,
  boundary estimates, objective value, and any numerical message.
- [ ] Calculate base-study power from the fitted distribution and the
  comparison's critical value.
- [ ] Scale noncentrality as `lambda(n) = lambda0 * n / n0` only after validating
  that `n0` is finite and positive.
- [ ] Remove PPE's incidental dependency on a fitted simulation-role schema:
  explicit comparisons between any two fitted model labels must work when
  `estimateSimulation = FALSE`.
- [ ] Make `method = "distribution_mle"` the documented default and retain the
  current algorithm as `method = "exceedance"`.
- [ ] Give the legacy method fields explicit names such as
  `threshold_exceedance_probability` and `threshold_implied_ncp`.
- [ ] Test the MLE against a direct likelihood grid and against fixed
  noncentral-chi-square simulation fixtures over weak, moderate, and strong
  effects.

**Acceptance:** Changing the plotting threshold changes calculated power but
does not change the distribution-MLE noncentrality estimate for the same
comparison data.

## Task 6: Add model-based PPE uncertainty

**Files:**

- Modify: `R/sse-ppe.R`
- Modify: `R/plot-sse.R`
- Modify: `tests/testthat/test-sse-ppe.R`

- [ ] Implement a parametric bootstrap that draws the observed number of
  positive test statistics from the fitted noncentral chi-square model and
  refits noncentrality.
- [ ] Give the bootstrap a deterministic keyed seed independent of simulation
  and estimation task scheduling.
- [ ] Report bootstrap failures and the number of successful bootstrap fits.
- [ ] Form percentile intervals for noncentrality and map each bootstrap draw
  through sample-size scaling to form power bands.
- [ ] Label these intervals `model_based` and keep them distinct from empirical
  binomial intervals.
- [ ] Permit `bootstrapSamples = 0` to skip the cost without changing the point
  estimate.
- [ ] Test reproducibility across sequential and supported parallel execution.

**Acceptance:** Plot data expose point estimate, confidence limits, bootstrap
replicate count, bootstrap failures, and interval type; no field is presented
as general model or effect-size uncertainty.

## Task 7: Add PPE distribution diagnostics

**Files:**

- Create: `R/plot-sse-ppe-diagnostics.R`
- Create: `tests/testthat/test-plot-sse-ppe-diagnostics.R`
- Modify: `NAMESPACE`
- Create: `_pkgdown.yml` (the package does not currently have one)

- [ ] Export `plotSSEPpeDiagnostics()`.
- [ ] Plot the empirical CDF of positive test statistics against the fitted
  noncentral chi-square CDF.
- [ ] Add a quantile-quantile panel with an identity line.
- [ ] Add pointwise parametric-bootstrap envelopes and label them as pointwise,
  not simultaneous.
- [ ] Calculate a prespecified discrepancy statistic, preferably
  Cramér--von Mises, and a parametric-bootstrap diagnostic p-value.
- [ ] Display the number and fraction of zero and negative statistics outside
  the fitted support.
- [ ] Return all diagnostic data in the `ggplot` object or an attached
  documented attribute so results are auditable without reading pixels.
- [ ] Treat diagnostics as evidence about approximation adequacy, not as a
  mechanical pass/fail certification.

**Acceptance:** Fixtures generated from the assumed distribution generally lie
inside the envelope, while deliberately misspecified mixtures visibly depart
and yield small bootstrap diagnostic p-values under a fixed test seed.

## Task 8: Validate sample-size proportionality

**Files:**

- Create: `R/sse-ppe-validation.R`
- Create: `tests/testthat/test-sse-ppe-validation.R`
- Modify: `NAMESPACE`
- Modify: `_pkgdown.yml`

- [ ] Export `validateSSEPpeScaling()` accepting completed SSE objects from at
  least two directly simulated study sizes and one explicit comparison.
- [ ] Verify that comparison definitions, parameter source, model labels, and
  design metadata are compatible across runs.
- [ ] Fit noncentrality independently at each study size.
- [ ] Report `lambda / n`, its uncertainty, and deviations from a through-origin
  linear relationship.
- [ ] Plot observed noncentralities with the proportional scaling line and
  uncertainty intervals.
- [ ] Require at least three study sizes before calculating a lack-of-fit
  statistic; with two sizes, provide only the descriptive ratio comparison.
- [ ] Do not pool raw replicates across study sizes.

**Acceptance:** The helper detects a deliberately nonlinear fixture and gives a
clear warning that single-size PPE extrapolation is unsupported over that
range.

## Task 9: Correct Monte Carlo uncertainty for parameter summaries

**Files:**

- Modify: `R/sse-helpers.R`
- Modify: `R/sse-methods.R`
- Modify: `tests/testthat/test-sse-statistical-contract.R`
- Create: `tests/testthat/test-sse-parameter-summary.R`

- [ ] Construct replicate-level absolute and relative errors before summary
  aggregation.
- [ ] Calculate bias MCSE from the sample SD of absolute errors divided by the
  square root of their effective count.
- [ ] Calculate relative-bias MCSE from the sample SD of replicate relative
  errors, including varying-truth runs after excluding zero truths.
- [ ] Ensure every standard error is nonnegative, including negative fixed
  truths.
- [ ] Report normal-approximation intervals only when the effective count and
  variance are sufficient; otherwise return `NA` with the effective count.
- [ ] Add bootstrap intervals as an optional future method, not as an unplanned
  dependency of this task.
- [ ] Replace ambiguous `rse` with `mcse_relative_bias`; preserve `rse` as a
  documented superseded output alias for one release if compatibility requires
  it. Because it is a table field rather than a callable argument, do not imply
  that lifecycle can warn when the field is read.
- [ ] Add `mcse_bias`, `ci_bias_lower`, `ci_bias_upper`,
  `mcse_relative_bias`, `ci_relative_bias_lower`, and
  `ci_relative_bias_upper` to tidy summaries.

**Acceptance:** Fixed positive truth, fixed negative truth, varying truth,
zero-truth exclusions, missing estimates, and one-observation cases all have
analytically checked expected results.

## Task 10: Add parameter-draw adequacy summaries

**Files:**

- Create: `R/sse-parameter-diagnostics.R`
- Create: `tests/testthat/test-sse-parameter-diagnostics.R`
- Modify: `R/run-sse.R`
- Modify: `R/sse-methods.R`
- Modify: `NAMESPACE`
- Modify: `_pkgdown.yml`

- [ ] Export `parameterDrawSummary()` and `plotSSEParameterDraws()`.
- [ ] Report realized mean, SD, median, central quantiles, minimum, maximum, and
  finite count for every generating parameter.
- [ ] For covariance mode, persist enough target metadata to compare realized
  THETA means/SDs and OMEGA means with the intended local approximation.
- [ ] For `joint`, report raw-scale OMEGA mean drift and empirical dependence
  for covered parameters without claiming exact covariance preservation.
- [ ] For `independent_iw`, report each block's binding degrees of freedom and
  intended versus realized diagonal dispersion.
- [ ] Check every generated OMEGA block for positive-definiteness.
- [ ] Compare THETA draws with finite model bounds where those bounds can be
  recovered reliably; report out-of-domain counts but do not silently truncate
  or resample.
- [ ] Store diagnostics in the returned object and `sse_summary.rds`.
- [ ] Warn only on prespecified actionable conditions; always retain the full
  table.

**Acceptance:** The diagnostics expose known mean inflation in a 1-by-1 `joint`
fixture, nonbinding over-dispersion in an `independent_iw` fixture, and an
out-of-bound THETA fixture.

## Task 11: Documentation, migration, and package validation

**Files:**

- Modify: `docs/sse-technical-reference.md`
- Modify: `README.md`
- Modify: `vignettes/runSSE.Rmd`
- Modify: relevant roxygen comments and generated `.Rd` files
- Modify: `_pkgdown.yml`
- Modify: `NEWS.md`

- [ ] Document the estimand for empirical power, the paired denominator, and
  all Monte Carlo uncertainty fields.
- [ ] Document starting-value policy as part of the numerical experiment.
- [ ] Document distribution-MLE PPE, negative-statistic handling, bootstrap
  uncertainty, diagnostics, and proportional-scaling assumptions.
- [ ] State that an explicit comparison does not prove nesting or justify its
  null distribution.
- [ ] Explain why covariance draw modes remain approximations and why no
  natural-coordinate compatibility sampler is being added.
- [ ] Add lifecycle documentation for `randomEstimationInits`; describe legacy
  implicit comparisons, `rse`, and `method = "exceedance"` as compatibility or
  superseded behavior without emitting misleading lifecycle warnings.
- [ ] Add NEWS bullets for every user-visible API or output change.
- [ ] Regenerate documentation with the package's pinned toolchain and inspect
  the diff for unrelated generated changes.
- [ ] Run targeted test files after each task, then `devtools::test()`,
  `devtools::check()`, and `pkgdown::check_pkgdown()`.
- [ ] Render the technical reference and vignette and verify formulas, tables,
  links, and plot captions.

**Acceptance:** Documentation identifies which uncertainty is empirical,
Monte Carlo, or model-based; all new exports appear in the reference index;
tests and package checks introduce no unexplained regressions.

---

## Statistical validation matrix

| Concern | Required validation |
| --- | --- |
| OFV direction | Hand-calculated full/reduced fixture with known sign |
| Paired denominator | Asymmetric missingness and filtering fixture |
| Empirical power CI | Exact binomial calculation from reported counts |
| Parameter bias MCSE | Hand calculation for fixed, negative, and varying truths |
| Initialization sensitivity | Same datasets and truths, different fit starts |
| PPE NCP MLE | Direct likelihood grid and simulated known-NCP fixtures |
| PPE bootstrap | Fixed-seed coverage experiment and failure accounting |
| PPE distribution diagnostic | Correct-model and misspecified-mixture fixtures |
| Sample-size scaling | Linear and deliberately nonlinear multi-size fixtures |
| `independent_iw` draws | Mean, binding/nonbinding SD, and positive-definiteness |
| `joint` draws | Positive-definiteness, mean drift, and empirical dependence |
| Resume/recompute | Statistical definitions and RNG provenance unchanged |

Monte Carlo validation tolerances must be selected before looking at results,
be wide enough for supported platforms, and be justified from the expected
Monte Carlo standard error rather than chosen to make a test pass.

## Risks and mitigations

### Distribution-MLE PPE can look authoritative when its assumptions fail

Mitigation: require explicit comparisons, report negative differences, provide
ECDF/QQ diagnostics, and support direct multi-size validation. Never describe a
bootstrap band as covering model misspecification.

### Excluding failed fits can bias operating characteristics

Mitigation: always expose attempted and paired-evaluable counts, summarize
failure mechanisms, and encourage sensitivity analyses that classify failed
fits according to a prespecified scientific rule rather than silently dropping
them.

### Starting at the generating truth can overstate numerical performance

Mitigation: keep stored model starts as defaults, record starts by role, and
make same-dataset sensitivity comparisons straightforward.

### Parameter-count degrees of freedom can be wrong

Mitigation: mark inferred degrees of freedom, warn when used for PPE, and make
explicit comparison objects the documented workflow. Boundary and constrained
hypotheses require user justification or empirical null calibration.

### Multiple comparisons can inflate the overall error rate

Mitigation: label unadjusted results as marginal, retain replicate-level test
statistics so prespecified joint decision rules can be evaluated, and do not
apply an implicit multiplicity correction.

### Parametric bootstrap is computationally costly

Mitigation: bootstrap only the inexpensive fitted reference distribution, use
keyed reproducible RNG, allow zero bootstrap samples, and persist fitted PPE
objects where appropriate.

### API migration can silently change historical results

Mitigation: preserve explicit legacy methods for at least one release, attach
method and comparison metadata to every output, test recomputation of old run
directories, and announce default changes in NEWS.

## Out of scope

- NONMEM control-stream, `$PRIOR`, NWPRI, TNPRI, MSF, table, and covariance-file
  support.
- Exact PsN random-number or output-file reproduction.
- Natural-scale multivariate-Normal OMEGA draws that can be indefinite.
- Automatic proof that models are nested or that their OFVs are comparable.
- Automatic selection of boundary-mixture null distributions.
- Treating covariance approximations as Bayesian posteriors.
- Silently clipping, truncating, repairing, or resampling scientifically invalid
  parameter draws.
- Changes to the underlying `nlmixr2est` fitting algorithms.

## Definition of done

The work is complete when:

1. every reported comparison has an auditable definition and paired
   denominator;
2. empirical probabilities and parameter biases include correctly named Monte
   Carlo uncertainty;
3. starting-value policy is role-specific, reproducible, and separable from the
   simulated-data distribution;
4. distribution-MLE PPE estimates one noncentrality parameter per comparison,
   with bootstrap uncertainty and distributional diagnostics;
5. multi-size SSE runs can assess the proportional-noncentrality assumption;
6. realized parameter draws have standard adequacy diagnostics;
7. legacy outputs remain reproducible through explicit compatibility options;
   and
8. documentation and tests state limitations without claiming that software
   can verify the scientific assumptions supplied by the user.

## References

1. Ueckert S, Karlsson MO, Hooker AC. Accelerating Monte Carlo power studies
   through parametric power estimation. *Journal of Pharmacokinetics and
   Pharmacodynamics*. 2016;43:223-234.
   [doi:10.1007/s10928-016-9468-y](https://doi.org/10.1007/s10928-016-9468-y).
2. Oehlert GW. A note on the delta method. *The American Statistician*.
   1992;46:27-29.
   [doi:10.1080/00031305.1992.10475842](https://doi.org/10.1080/00031305.1992.10475842).
3. Pinheiro JC, Bates DM. Unconstrained parametrizations for
   variance-covariance matrices. *Statistics and Computing*. 1996;6:289-296.
   [doi:10.1007/BF00140873](https://doi.org/10.1007/BF00140873).
4. Self SG, Liang K-Y. Asymptotic properties of maximum likelihood estimators
   and likelihood ratio tests under nonstandard conditions. *Journal of the
   American Statistical Association*. 1987;82:605-610.
   [doi:10.1080/01621459.1987.10478472](https://doi.org/10.1080/01621459.1987.10478472).
5. Barnard J, McCulloch R, Meng X-L. Modeling covariance matrices in terms of
   standard deviations and correlations, with application to shrinkage.
   *Statistica Sinica*. 2000;10:1281-1311.
   [Article](https://www3.stat.sinica.edu.tw/statistica/j10n4/10-4.htm).
6. Dosne A-G, Bergstrand M, Karlsson MO. An automated sampling importance
   resampling procedure for estimating parameter uncertainty. *Journal of
   Pharmacokinetics and Pharmacodynamics*. 2017;44:509-520.
   [doi:10.1007/s10928-017-9542-0](https://doi.org/10.1007/s10928-017-9542-0).
