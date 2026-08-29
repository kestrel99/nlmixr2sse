# Critical Review: SSE Parameter Uncertainty Design and Plan

**Date:** 2026-08-28  
**Artifacts reviewed:**

- [`2026-08-28-sse-parameter-uncertainty-design.md`](./2026-08-28-sse-parameter-uncertainty-design.md)
- [`2026-08-28-sse-parameter-uncertainty.md`](../plans/2026-08-28-sse-parameter-uncertainty.md)
- Relevant `nlmixr2est`, `nlmixr2sir`, `rxode2`, PsN, and `nlmixr2sse` source and tests

## Summary

**Verdict: Request changes.**

The overall direction is valuable: `parameterSource = "covariance"` should be able to propagate OMEGA uncertainty rather than varying theta parameters alone. The proposed log-Cholesky joint construction also addresses the important requirement that simulated OMEGA matrices remain positive definite.

The current design and implementation plan nevertheless have three blocking problems:

1. The proposed `"nwpri"` distribution is not NONMEM/PsN NWPRI.
2. Fixed or partially covered OMEGA blocks are not handled consistently or safely.
3. Common theta-only covariance cases can silently stop drawing theta or fail.

The statistical justification for making `"joint"` the default is also not yet strong enough. Its transformation bias can be substantial in precisely the poorly identified settings where parameter-uncertainty propagation is most important. Neither mode should be made the default until the semantics, fallbacks, and validation are corrected.

## Critical Issues (Blocking)

### 1. The proposed `"nwpri"` mode is not NONMEM/PsN NWPRI

The design defines `"nwpri"` as an independent, mean-centred textbook inverse-Wishart construction with

\[
\Psi=(\nu-p-1)\widehat\Omega.
\]

That is not the NONMEM NWPRI parameterization.

`rxode2` explicitly distinguishes NWPRI from a general inverse-Wishart density:

- `../rxode2/R/priorDensity.R`, lines 22-39
- `../rxode2/inst/include/rxode2prior.h`, lines 20-27 and 67-69
- `../rxode2/src/priorDensity.cpp`, lines 289-305
- `../rxode2/tests/testthat/test-prior-density.R`, lines 539-600

The tests specifically verify that the NWPRI and general inverse-Wishart densities differ. `rxode2` also already implements PRIOR-style simulation:

- `../rxode2/R/prior-sim.R`, lines 282-286, 388-420, and 593-608
- `../rxode2/src/cvPost.cpp`, lines 77-117

For example, its scalar construction is proportional to `nu * omega / gamma(nu / 2, scale = 2)`, and the matrix draw is likewise scaled according to NONMEM-style semantics. This is not the design's mean-centred `Psi = (nu - p - 1) * Omega` construction.

PsN does not infer OMEGAPD values from standard errors. `update_inits` asks users to supply blockwise degrees of freedom through `add_prior`:

- `../PsN/doc/update_inits_userguide.tex`, line 52
- `../PsN/bin/update_inits`, lines 139-148
- `../PsN/lib/model/problem.pm`, lines 306-316 and 449-452

PsN collects `omega_variance` values but does not use them. That is insufficient evidence for the design's assertion that PsN intended the proposed moment-matching construction. Describing the new procedure as completing an unfinished PsN intent is speculation.

Required resolution:

- Rename the proposed distribution to something such as `"independent_iw"` or `"iw_mean"`; or
- Implement genuine NWPRI semantics, preferably by using or validating against `rxode2::cvPost`, and expose user-supplied blockwise OMEGAPD/\(\rho\) values.

Both distributions could be offered, but they are distinct statistical procedures. Claims of PsN parity must be reserved for the genuine NWPRI implementation.

### 2. Fixed and partially covered OMEGA blocks are unsafe

`nlmixr2est` deliberately removes fixed OMEGA elements from `fit$cov`:

- `../nlmixr2est/R/foceiCov.R`, lines 35-37

The declared OMEGA block topology still contains those elements. The plan's joint implementation constructs `joint_names` from every entry in a selected OMEGA block and requires every name to exist in `fit$cov`. A correlated block containing one fixed component will therefore abort.

The planned independent implementation has the opposite error: once any diagonal in a block has a usable standard error, `.drawOmega()` draws the entire block. This can mutate fixed or uncovered diagonal and covariance elements, contradicting the design's promise to remove and restore fixed components.

Restoring fixed entries after drawing the rest is not generally safe. Replacing entries of a positive-definite draw with fitted cross-covariances can make the final matrix non-positive-definite.

For the first implementation, the safe policy is:

> If any element of a correlated OMEGA block is fixed or lacks covariance coverage, hold the entire block fixed and report why.

Drawing the estimable subspace of a partially fixed block would require a constrained parameterization that guarantees positive definiteness. It should be treated as a later extension, not implemented through element replacement.

### 3. Theta-only covariance is a supported nlmixr2 case

The design overstates the general availability of full theta-OMEGA covariance. `nlmixr2est` does use the expected natural-scale OMEGA names:

- `../nlmixr2est/R/foceiCov.R`, lines 77-84
- `../nlmixr2est/R/foceiCovAnalytic.R`, lines 708 and 771

Full covariance is nevertheless conditional:

- Users can request theta-only covariance using `covFull = FALSE`: `../nlmixr2est/R/foceiControl.R`, line 104.
- Full covariance is unavailable for some bounded transformations or when no theta rows exist: `../nlmixr2est/R/foceiCovFdFull.R`, lines 17-22.
- The installer retains the native theta-only covariance when required components are missing, non-finite, or non-positive-definite: `../nlmixr2est/R/foceiCovFdFull.R`, lines 37-40 and 69-80.

The plan creates a joint specification only when at least one OMEGA block is drawable. With a valid theta-only `fit$cov`, the joint specification is `NULL`, and the proposed resolution path returns fitted theta values rather than drawing them.

Required behavior:

- Always draw covered theta parameters.
- Draw only fully supported OMEGA blocks.
- Record uncovered OMEGA blocks as fixed.
- Treat partial covariance coverage as a supported degraded mode rather than an all-or-nothing failure.

There is also an OMEGA-only indexing bug in the planned `.drawJoint()` implementation. In R, `drawn[-seq_len(spec$nTheta)]` returns an empty vector when `nTheta == 0`, so an OMEGA-only joint draw is lost.

## Required Changes

### 4. Joint mode incorporates covariance only approximately

The claim that joint mode “preserves correlation” is too strong. The log-Cholesky transformation and delta method incorporate the fitted covariance through a first-order local approximation. After nonlinear back-transformation, they do not exactly preserve raw-scale means, variances, Pearson correlations, or the fitted asymptotic distribution.

The documentation should say that joint mode:

- Incorporates theta-OMEGA and OMEGA-OMEGA cross-covariance through a first-order delta approximation.
- Guarantees positive-definite OMEGA draws.
- Does not exactly reproduce raw-scale fitted moments.

Similarly, `fit$cov` is an estimated covariance matrix used in a local asymptotic normal approximation. It is not itself “the asymptotic sampling distribution.”

### 5. Joint-mode mean bias is not necessarily slight

For a scalar OMEGA, let \(\phi=\frac12\log\Omega\). The proposed centred Gaussian draw gives

\[
\frac{E(\Omega^*)}{\widehat\Omega}
=
\exp\left(\frac{\operatorname{SE}(\widehat\Omega)^2}
{2\widehat\Omega^2}\right).
\]

This implies approximately:

| OMEGA relative SE | Mean inflation |
| ---: | ---: |
| 20% | 2.0% |
| 50% | 13.3% |
| 100% | 64.9% |

The single 20% example in the design is not enough to characterize the bias as slight. Sparse or poorly identified studies can have much larger relative standard errors, and those are precisely the cases in which parameter uncertainty matters most.

If joint mode is to be the default, mean correction or calibration should not remain out of scope. At minimum, the mode should initially be opt-in with explicit diagnostics for the expected transformation bias.

### 6. Evidence for making joint mode the default is insufficient

A maximum correlation of 0.442 in one thinned `theo` example demonstrates that nonzero covariance can occur. It does not show that ignoring it materially distorts SSE operating characteristics.

Choosing a default should be based on downstream quantities, including:

- Type-I error and power.
- The distribution of likelihood-ratio statistics.
- Prediction and exposure distributions.
- Bias relative to bootstrap or raw-results sampling.
- Behavior across multiple model dimensions, designs, and uncertainty levels.

`nlmixr2sir` supports the premise that cross-covariance can matter because its proposal uses full `fit$cov` (`../nlmixr2sir/R/sir.R`, lines 82-111). It does not validate the proposed log-Cholesky Gaussian distribution: SIR samples on the natural scale, caps correlations by default, and rejects invalid OMEGA draws.

### 7. The independent inverse-Wishart claims need narrower language

Under the proposed textbook parameterization, the generated OMEGA draws are mean-centred at `omegaHat`. Calling them “exactly unbiased” is misleading: unbiasedness normally describes an estimator over repeated datasets, not merely the expectation of a generated distribution.

Similarly, choosing the minimum blockwise \(\nu\) only ensures that the constructed marginal diagonal standard deviations are not smaller than the requested diagonal values. It is not guaranteed to be conservative for nonlinear predictions, correlations, likelihood-ratio statistics, Type-I error, or power.

One \(\nu\) per block also:

- Ignores off-diagonal OMEGA standard errors.
- Ignores covariance among OMEGA estimators.
- Allows the least precise diagonal to inflate the entire block.

Use “mean-centred” and “matches or exceeds the selected diagonal marginal uncertainty” rather than “unbiased” and “conservative.”

### 8. The numerical Jacobian needs scale-aware perturbations

The planned absolute step `h = 1e-6` is scale-sensitive. It is too large for very small OMEGA components, too small for large components, and can perturb a near-boundary covariance matrix outside the positive-definite cone.

Use relative/adaptive steps with step-halving and a one-sided fallback, or implement analytic/automatic derivatives. Tests should cover:

- OMEGA scales from approximately `1e-8` to `1e2`.
- Correlations near their admissible boundaries.
- Ill-conditioned but positive-definite blocks.
- One-, two-, and three-dimensional blocks.

### 9. The planned degrees-of-freedom error is effectively unreachable

The proposed formula constructs `nu` as `p + 3 + positive_term`. Consequently, checking `nu <= p + 3` does not detect ordinary weak identification; it is reachable mainly through numerical underflow or invalid inputs.

If weak identification should produce a warning or error, define an explicit policy based on relative SE, the margin `nu - (p + 3)`, or a user-configurable threshold. Also validate the full fitted OMEGA block as positive definite before calling `rWishart`, with an error identifying the block and ETA names.

### 10. Theta constraints remain unaddressed

The current covariance mode already draws theta values from an unconstrained multivariate normal. Extending the mode to “full parameter uncertainty” makes this limitation more visible: bounded theta parameters and residual-error standard deviations can be drawn outside their valid domains.

This need not block the initial OMEGA work, but the documentation must state the limitation. A later design should consider transformed or truncated theta draws.

### 11. The implementation plan contradicts itself

The plan opens by mandating independent NWPRI-only behavior and later says not to redesign that decision. Subsequent tasks instead add both modes and make joint the default. These are mutually exclusive implementation instructions.

The plan introduction and its “design settled” section must be replaced with the final agreed architecture before implementation begins.

### 12. Existing tests and documentation will become false

The plan does not update all existing assertions and user-facing statements that OMEGA remains fixed:

- `tests/testthat/test-run-sse.R`, lines 168-199
- `tests/testthat/test-run-sse.R`, lines 878-888
- `R/run-sse.R`, lines 289-292
- `vignettes/runSSE.Rmd`, lines 200-220
- `README.md`, line 97
- `man/runSSEControl.Rd`, lines 34-35
- `NEWS.md`, line 29

The claim that the full suite will retain the same expected failures is therefore not credible. Once OMEGA drawing works, the current integration assertion that OMEGA has one unique value must fail.

### 13. `covarianceDraw` belongs to the reproducibility contract

Changing from joint to independent draws changes all simulation parameters generated from the same seed. The selected covariance-draw mode must be persisted in `run_info` and checked during resume and `addModels` operations.

Current validation checks `parameterSource`, `randomEstimationInits`, and `updateFix`, but not the proposed mode:

- `R/run-sse.R`, lines 178-198
- `R/recompute-sse.R`, lines 109-135

Backward compatibility also needs a policy for completed runs created before `covarianceDraw` existed.

## Relative Assessment of the Modes

| Criterion | Joint delta mode | Proposed independent IW | True NWPRI |
| --- | --- | --- | --- |
| Uses theta-OMEGA covariance | Yes, approximately | No | No |
| Uses between-block covariance | Yes, approximately | No | No |
| Produces positive-definite OMEGA | Yes | Yes | Yes |
| Preserves raw-scale OMEGA mean | No, unless corrected | Yes by construction | Not under the proposed mean formula |
| Needs complete full covariance | Needs an explicit fallback | Less dependent on it | No, if prior parameters are supplied |
| PsN/NONMEM parity | No | No | Yes |
| Primary strength | Retains estimated dependence | Simple, positive-definite, mean-centred blocks | Compatibility with NONMEM/PsN prior workflows |
| Primary weakness | Transformation approximation and bias | Discards dependence and reduces a block to one df | Requires meaningful prior df and does not represent full sampling covariance |
| Best role | Full covariance is reliable and dependence matters | Sensitivity analysis using independent mean-centred blocks | Reproducing NONMEM/PsN prior simulation semantics |

The most defensible API would distinguish the procedures explicitly:

- `"joint_delta"`: transformed joint asymptotic approximation.
- `"independent_iw"`: the proposed mean-centred independent inverse-Wishart construction.
- `"nwpri"`: genuine NONMEM-compatible NWPRI behavior.

## Strong Suggestions

### Distributional validation

Add oracle-based tests rather than testing implementation internals alone:

- Validate genuine NWPRI draws against `rxode2::cvPost` moments or distributions.
- Compare joint-delta results with parametric bootstrap or canonical raw-results draws.
- Assess downstream SSE Type-I error and power, not only parameter moments.

### Required edge-case tests

Add explicit tests for:

- Theta-only `fit$cov`.
- OMEGA-only covariance with zero theta parameters.
- A correlated OMEGA block containing a fixed ETA variance.
- Missing covariance entries for an otherwise declared block.
- One-, two-, and three-dimensional block ordering.
- Near-boundary and ill-conditioned OMEGA matrices.
- Accurate drawn/fixed parameter partition reporting.
- Resume and `addModels` requests with a changed draw mode.

## Recommended Revision Sequence

1. Decide whether `"nwpri"` means actual NONMEM NWPRI or rename the proposed construction.
2. Define explicit theta-only and partial-coverage fallbacks.
3. Hold partially fixed correlated blocks fixed in the initial implementation.
4. Add mean calibration or opt-in status for joint draws.
5. Replace “preserves correlation,” “exactly unbiased,” and “conservative” with narrower claims.
6. Resolve the contradictory opening of the implementation plan.
7. Expand tests, documentation, and resume compatibility before selecting a default.

Once these issues are addressed, joint-delta drawing is likely the more informative general-purpose approximation when a trustworthy full covariance is available. A true NWPRI mode remains useful for NONMEM/PsN compatibility, while the proposed mean-centred independent inverse-Wishart mode is better presented as a distinct sensitivity-analysis option.
