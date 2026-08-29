# SSE Parameter Uncertainty Design

## Goal

Make `runSSE(control = runSSEControl(parameterSource = "covariance"))` simulate
with uncertainty in **both** fixed effects and random effects, using the
covariance that current `nlmixr2est` reports over THETA and OMEGA jointly. This
restores a currently broken code path.

Two draw modes are provided, selected by `runSSEControl(covarianceDraw = )`:

| Mode | Description |
| --- | --- |
| `"independent_iw"` (**default**) | THETA multivariate-Normal and OMEGA mean-centred inverse-Wishart per block, drawn **independently**, discarding the cross-terms. No transformation bias. |
| `"joint"` (opt-in) | THETA and OMEGA drawn **together** from the full covariance on a log-Cholesky-transformed scale, incorporating the estimated THETA↔OMEGA covariance through a first-order delta approximation. Positive-definite by construction, but inflates OMEGA means when relative standard errors are large. |

**Neither mode is NONMEM `$PRIOR NWPRI`.** An earlier draft of this design
called the independent mode `"nwpri"` and claimed PsN parity for it. That was
wrong: `rxode2` documents that NONMEM's NWPRI does *not* reduce to a textbook
inverse-Wishart — its degrees-of-freedom convention is `d_W = rho + n + 1` with
scale `rho * Psi`, giving a different density and gradient
(`rxode2/R/priorDensity.R`). PsN also never derives `OMEGAPD` from standard
errors; `update_inits -add_prior` requires the user to supply blockwise degrees
of freedom. No parity claim is made for either mode here. A genuine NWPRI mode
remains possible later, and would be built on `rxode2`'s existing
implementation rather than re-derived.

`"independent_iw"` is the default for now because it carries no transformation
bias. `"joint"` is opt-in, and warns when the bias is expected to be material.
Whether `"joint"` should eventually become the default is deliberately left
open — see *Choosing a default* below.

## Background

### The immediate defect

`parameterSource = "covariance"` is broken against current `nlmixr2est`.

`nlmixr2est` now **can** report a joint fixed/random-effect covariance in
`fit$cov` — for FOCEI `"r"`, `"s"`, `"r,s"`, `"analytic"` and SAEM `"fim"`,
`"sa"` — whenever the full-covariance calculation and installation succeed. It
is **not** guaranteed: `foceiControl(covFull = FALSE)` requests theta-only
directly, some bounded transformations cannot produce it, and the installer
falls back to the native theta-only covariance when components are missing,
non-finite, or non-positive-definite. See *Coverage policy* below, which
treats theta-only as a supported configuration rather than a failure.

When it is reported, the matrix carries OMEGA rows named `om.<eta>` (variance)
and `cov.<eta_i>.<eta_j>` (covariance) alongside the theta rows:

```
dimnames(fit$cov)[[1]]
#> "tka" "tcl" "tv" "add.sd" "om.eta.ka"
```

`.alignedCovariance()` assumes `fit$cov` is theta-only and aborts on the first
name absent from `fit$theta`:

```
Error in .alignedCovariance(fit) :
  `fit$cov` contains theta name "om.eta.ka" that are not present in `fit$theta`.
```

This reproduces with the exact model and control settings in the package
README, so covariance-mode SSE fails for ordinary usage. The package test suite
does not catch it: the one integration test covering covariance mode is masked
by an unrelated `cli` defect in `nlmixr2utils:::.abortRawResults()` (an
`.envir` problem that turns a real error into `skip()`).

### The opportunity

The same `fit$cov` OMEGA rows that break the current code are exactly the
information needed to draw OMEGA with uncertainty. Until recently that
information was not available, which is why the mode was documented as
theta-only ("OMEGA and SIGMA stay fixed"). It is available now.

## Relationship to PsN

PsN's `sse` tool does not construct priors itself; it detects a `$PRIOR NWPRI`
or `TNPRI` record already present in the input model and lets NONMEM sample
from it. The prior is built beforehand by `update_inits -add_prior=df1,df2,...`,
implemented in `add_prior_distribution()` (PsN `lib/model/problem.pm:306`),
which reads a previous run's estimates and covariance matrix from output and
generates:

| Record | Source |
| --- | --- |
| `THETAP` | final THETA estimates (FIX) |
| `THETAPV` | leading `ntheta × ntheta` block of the `.cov` matrix |
| `OMEGAP` | final OMEGA estimates, preserving BLOCK structure |
| `OMEGAPD` | user-supplied degrees of freedom, one integer per OMEGA block |

NWPRI's statistical structure is a product of independent priors —
`p(theta, Omega) = p(theta) · p(Omega)` — with THETA multivariate Normal and
OMEGA inverse-Wishart. There is no cross-covariance term between the blocks.

`"independent_iw"` shares that factorization but not that density. How the two
modes compare:

| | PsN / NWPRI | This design |
| --- | --- | --- |
| Prior construction | automated, separate `update_inits` step | automatic, inside `runSSE()` |
| THETA prior variance | theta block of `$COV` | same |
| OMEGA centre | final OMEGA estimates | same (`fit$omega`) |
| OMEGA density | NONMEM `d_W = rho + n + 1`, scale `rho * Psi` | **textbook mean-centred inverse-Wishart** — a different distribution |
| Independence of blocks | yes | `"independent_iw"`: same. `"joint"` (opt-in): THETA↔OMEGA covariance incorporated to first order |
| OMEGA df | user supplies, one per block | `"independent_iw"`: **moment-matched per block from OMEGA SEs**. Not applicable to `"joint"` |
| SIGMA block | none (`neps=0` hardcoded) | none (residual error is theta) |
| Fixed parameters | "known to not work when one or more parameters are fixed" | designed for from the start |

**Neither mode implements NWPRI, and neither claims PsN parity.** The table
above describes what PsN does; the differences below are departures, not
equivalences.

**1. `"independent_iw"` shares NWPRI's independence but not its density.**
NWPRI factorizes as `p(θ)·p(Ω)`, and `"independent_iw"` does too, but the
OMEGA factor is a textbook mean-centred inverse-Wishart, not NONMEM's
`d_W = rho + n + 1` / `rho * Psi` construction. Same factorization, different
distribution.

**2. `"independent_iw"`'s degrees of freedom are derived, not user-supplied.**
PsN requires the user to supply one integer per block via
`update_inits -add_prior`. Deriving `nu` from the reported OMEGA standard
errors is a convenience this design adds, not a completion of PsN's intent.
(PsN does collect `@omega_variance` from the `$COV` diagonal at OMEGA
positions and never use it — `problem.pm:385-388` — but dead code is not
evidence of intent, and an earlier draft of this design over-read it.)

**3. `"joint"` uses the cross-terms both PsN and `"independent_iw"` discard.**
NWPRI's independence is a conjugacy convenience for specifying priors, not a
claim that estimates are uncorrelated. Since `fit$cov` reports the covariance,
and the log-Cholesky transform removes the positive-definiteness objection to
using it, `"joint"` incorporates it — approximately, to first order. Supporting
context: SIR (Dosne et al. 2016) uses the full variance–covariance matrix as
its proposal with correlations intact, and notes that mis-specified
correlations are hard to recover from. Note SIR samples on the natural scale
and rejects invalid draws, so it validates the *premise* that cross-covariance
matters, not this design's particular transformed-Gaussian construction.

The inverse-Wishart's known limitations bear on `"independent_iw"`: a single
degrees-of-freedom parameter governs the uncertainty of every variance in a
block, ignoring off-diagonal standard errors and covariance among the OMEGA
estimators, and it induces an a priori dependence between variances and
correlations (Barnard, McCulloch & Meng 2000; Alvarez, Niemi & Simpson 2014).
The single-`nu` consequence is quantified under Statistical Design.

**On SIGMA.** nlmixr2 has no separate SIGMA matrix. Residual-error parameters
(`add.sd`, `prop.sd`, …) are ordinary thetas and carry full standard errors in
both `fit$cov` and `parFixedDf`. They are therefore drawn by the THETA step
with no special handling, and residual-error uncertainty is fully represented.
`fit$sigma` returning `named numeric(0)` is a legacy accessor artifact
(`.sigma()` still reads the deprecated `x$uif`), not missing uncertainty. PsN's
automated path likewise emits no `SIGMAP` block.

## Chosen Approach

Extend the existing draw machinery in `.resolveCovarianceParameterSets()` and
`.alignedCovariance()`. Both modes share the same partitioning of `fit$cov`
into THETA and OMEGA entries — the step that fixes the reported defect — and
differ only in how the draw is taken.

**`covarianceDraw = "joint"` (opt-in).** Transform each OMEGA block to an
unconstrained log-Cholesky vector, carry `fit$cov` onto that scale with the
delta method, and draw THETA and the transformed OMEGA **together** from one
multivariate Normal. Back-transforming reconstructs OMEGA. Because the
back-transform exponentiates the Cholesky diagonal, every draw is
positive-definite by construction — no rejection sampling.

**`covarianceDraw = "independent_iw"`.** THETA multivariate Normal from the
theta sub-block, OMEGA inverse-Wishart per block mean-centred on `fit$omega`
with degrees of freedom moment-matched to the reported OMEGA standard errors,
the two drawn independently. The OMEGA draw itself is delegated to
`rxode2::cvPost()` rather than re-implemented.

Backward compatibility is explicitly **not** a goal: OMEGA is drawn whenever
the covariance supports it. Same-seed draws will differ from previous releases.

### Choosing a default

**The default is `"independent_iw"`, and that choice is provisional.** It is
provisional in a specific, bounded sense: the *implementation* has a settled
default that tests, documentation, and `run_info` must all agree on, while the
*question* of whether `"joint"` should eventually take its place stays open
pending the evidence listed below.

`"independent_iw"` was chosen because it carries no transformation bias. The
case for `"joint"` is suggestive but not yet sufficient, and its bias is
largest in exactly the poorly-identified settings where uncertainty propagation
matters most — so it would be the wrong thing to apply silently to users who
never chose it.

What argues for `"joint"`:

- A naive joint Normal draw over raw variance elements would not be
  positive-definite, which is what originally motivated the independent
  structure. The log-Cholesky transform removes that objection, so
  positive-definiteness is no longer a reason to discard the cross-terms.
- `fit$cov` is an estimated covariance used in a local asymptotic-normal
  approximation. Zeroing its THETA↔OMEGA block discards information the fit
  actually reports, and for any quantity depending on both (a CV, a
  concentration prediction, ΔOFV itself) the delta-method variance loses its
  `2∇g_θ′Σ_θΩ∇g_Ω` term.
- **Measured on this package's own reference model**, max |corr(θ, ω)| is
  0.074 for the rich `theo_sd` design but **0.442** when thinned to three
  samples per subject. The population mean and the between-subject variance of
  the same parameter compete to explain sparse data, so their estimators
  correlate strongly. Sample-size studies are that regime.
- `parameterSource = "rawres"` draws bootstrap parameter vectors, which carry
  the joint dependence automatically. An independent-block `"covariance"` mode
  represents uncertainty differently from `"rawres"` for the same fit.

Why that is **not yet enough** to justify a default:

- One maximum correlation in one thinned example shows that non-zero
  covariance occurs. It does not show that ignoring it materially distorts SSE
  operating characteristics.
- `"joint"`'s mean inflation grows as `exp(SE²/2Ω̂²)` — 2.0% at 20% relative
  SE, 13.3% at 50%, **64.9% at 100%** (derivation under Risks). Poorly
  identified variance components routinely sit at the high end.

A defensible default requires evidence on downstream quantities, not parameter
moments alone: Type-I error and power, the ΔOFV distribution, exposure
predictions, and bias against bootstrap or raw-results sampling, across
several model dimensions, designs, and uncertainty levels. Until that exists,
the honest position is that `"joint"` is the more informative approximation
*when the full covariance is trustworthy and uncertainty is moderate*, and
`"independent_iw"` is a simpler, positive-definite, mean-centred sensitivity
analysis.

## Statistical Design — `"joint"`

### Log-Cholesky transform

Each OMEGA block is mapped to an unconstrained vector:

- **Forward** (Ω → φ): `L = t(chol(Ω))`; take `L`'s lower-triangular elements
  in column-major order, replacing the diagonal entries with their logs.
- **Inverse** (φ → Ω): rebuild lower-triangular `L`, exponentiating the
  diagonal, then `Ω = L Lᵀ`.

Exponentiating the diagonal makes `L`'s diagonal strictly positive, so `L Lᵀ`
is positive-definite for **every** value of φ. This is what removes the need
for rejection sampling.

### Carrying the covariance onto the transformed scale

Let `ω` be the free OMEGA elements on the natural scale and `φ = g(ω)` the
log-Cholesky vector. With `J = ∂φ/∂ω`, the delta method gives:

```
Σ_T = B Σ Bᵀ,    B = blockdiag(I, J)
```

so `Var(φ̂) ≈ J Σ_ωω Jᵀ` and, crucially, `Cov(θ̂, φ̂) ≈ Σ_θω Jᵀ` — the
cross-block is carried through rather than discarded.

`J` is computed by central finite differences over the forward transform. The
number of free OMEGA elements is small (typically well under ten), so this is
cheap, and no new package dependency is required.

**The step must be scale-aware.** A fixed absolute step (`h = 1e-6`) is wrong
for OMEGA: too large for components near `1e-8`, too small for components near
`1e2`, and — worst — capable of perturbing a near-boundary covariance matrix
outside the positive-definite cone, where `chol()` fails. Use a relative step,
`h_j = max(|x_j|, x_typical) * eps^(1/3)`, with:

- **step halving** on failure: if either perturbed point cannot be transformed
  (a `chol()` failure), halve the step and retry, to a bounded number of
  attempts;
- **one-sided fallback**: if the step cannot be made to work in both
  directions, fall back to a forward difference on the side that succeeds;
- **abort with the block and eta names** if neither direction works, rather
  than returning a silently wrong Jacobian.

### Drawing

One `chol(Σ_T)` and one `rnorm()` per replicate over the stacked
`(θ, φ)` vector, then back-transform each OMEGA block. Blocks are transformed
independently of one another (off-block OMEGA elements are structurally zero),
but all of them are stacked with THETA into a **single** joint draw, so both
THETA↔OMEGA and within-block OMEGA↔OMEGA covariances are incorporated to first order.

### Verified numerically

Two thetas plus a 2×2 OMEGA block, with `corr(tka, om₁₁) = 0.45` and
`corr(add, om₂₂) = −0.35` deliberately imposed; 40,000 draws:

| Quantity | Recovered | Target |
| --- | --- | --- |
| Positive-definite draws | 40000 / 40000 | all |
| `corr(tka, om₁₁)` | 0.445 | 0.45 |
| `corr(add, om₂₂)` | −0.342 | −0.35 |
| `SD(om₁₁)`, `SD(om₂₂)` | 0.0616, 0.0321 | 0.060, 0.030 |
| `E[om₁₁]`, `E[om₂₂]` | 0.3065, 0.1256 | 0.30, 0.12 |

Correlations and standard deviations are recovered closely. The means come out
about 2% high — see Risks.

## Statistical Design — `"independent_iw"`

### Inverse-Wishart parameterization

For `Omega ~ InvWishart(Psi, nu)` over a `p × p` block:

```
E[Omega_ii]   = Psi_ii / (nu - p - 1)
Var(Omega_ii) = 2 * Psi_ii^2 / ((nu - p - 1)^2 * (nu - p - 3))
```

Setting `Psi = (nu - p - 1) * Omega0` makes `E[Omega] = Omega0` exactly and
collapses the variance to:

```
Var(Omega_ii) = 2 * Omega0_ii^2 / (nu - p - 3)
```

Solving against the reported standard error
`SE_i = sqrt(fit$cov["om.<eta_i>", "om.<eta_i>"])` gives a closed form:

```
nu_i = p + 3 + 2 * (Omega0_ii / SE_i)^2
```

Take `nu = min(nu_i)` across the block's diagonal elements.

**Verified numerically.** For a 2×2 block with `Omega0 = [[0.30, 0.05], [0.05,
0.12]]` and target SEs `(0.06, 0.03)`, 200,000 draws gave empirical
`E[Omega_ii] = (0.300, 0.120)` against target `(0.300, 0.120)`, and empirical
`SD = (0.0751, 0.0300)` against closed-form `(0.0750, 0.0300)`.

**Known limitation, by construction.** The inverse-Wishart has a single `nu`
per block, so only the binding (minimum-`nu`) element matches its reported SE
exactly. In the verification above the non-binding element received `SD =
0.075` against its `0.06` target — 25% wider. Taking `min(nu_i)` makes every
non-binding element *more* dispersed than its reported SE, never less, so the
error is wider for the constructed diagonal marginal, never narrower. That is
not the same as being conservative for anything derived from the draws. This
must be documented in the user
guidance.

### Sampling mechanics

Delegated to `rxode2::cvPost(nu, omega * (nu - p - 1)/nu)`. The pre-scaling cancels `cvPost`'s `Psi = nu * omega` convention, giving a mean-centred draw. Positive-definiteness is guaranteed by construction, so no rejection sampling is required.

### Block structure

OMEGA blocks are derived from `fit$iniDf` as connected components over the
off-diagonal rows (`neta1 != neta2`). Diagonal-only etas form 1×1 blocks. Each
block gets its own `nu` and its own independent draw, matching PsN's
per-block `OMEGAPD` granularity.

## Coverage policy (shared by both modes)

`fit$cov` does not always cover every parameter, and the two cases must be
handled independently of each other.

### Partial covariance coverage is a supported degraded mode

Coverage is decided **per parameter and per OMEGA block**, never all-or-nothing:

- **Always draw covered thetas.** A theta present in `fit$cov` is drawn even
  when no OMEGA block is drawable.
- **Draw only fully covered OMEGA blocks.**
- **Record everything else as fixed** in `parameterPartition$fixed`, with the
  reason available in the run record.

A theta-only `fit$cov` is a *supported configuration*, not a failure:
`foceiControl(covFull = FALSE)` requests it directly
(`nlmixr2est/R/foceiControl.R`), full covariance is unavailable for some
bounded transformations, and the installer falls back to the native theta-only
covariance whenever the full-covariance components are missing, non-finite, or
non-positive-definite (`nlmixr2est/R/foceiCovFdFull.R`). In that case the run
draws thetas and holds all of OMEGA fixed — which is exactly the mode's
previous documented behaviour.

### A block with any fixed or uncovered element is held entirely fixed

`nlmixr2est` deliberately drops fixed OMEGA elements from `fit$cov`
(`.foceiOmegaPairs()` in `nlmixr2est/R/foceiCov.R` filters on `.omegaFixed()`),
while the declared block topology still contains them. So a correlated block
with one fixed component has entries missing from `fit$cov`.

**Policy: if any element of a correlated OMEGA block is fixed or lacks
covariance coverage, hold the whole block at its fitted values and record why.**

Two rejected alternatives, and why:

- *Draw the block anyway once any diagonal has an SE* — mutates fixed and
  uncovered elements, contradicting the guarantee that fixed parameters stay
  fixed.
- *Draw the free sub-matrix and restore fixed entries afterwards* — splicing
  fitted cross-covariances back into a positive-definite draw can produce a
  non-positive-definite matrix. Element replacement is not a safe repair.

Drawing the estimable subspace of a partially fixed block needs a constrained
parameterization that guarantees positive-definiteness. That is a later
extension, out of scope here.

A 1×1 block whose single element is fixed or uncovered is simply not drawn —
the same rule, with no special case.

### Reporting

`runInfo$parameterPartition$drawn` / `$fixed` is extended to cover OMEGA
entries, so the run record states exactly what varied and what did not.

## Component Design

| Function | Change |
| --- | --- |
| `.alignedCovariance()` | Partition `fit$cov` dimnames into theta and OMEGA subsets instead of aborting on OMEGA names. Return both blocks plus the fixed/uncovered sets. Shared by both modes. |
| `.omegaBlocks()` *(new)* | Derive OMEGA block structure as connected components over the off-diagonal entries. Shared by both modes. |
| `.omegaToPhi()` / `.phiToOmega()` *(new)* | Log-Cholesky forward and inverse transforms for one OMEGA block. |
| `.numericJacobian()` *(new)* | Central-difference Jacobian of a vector-valued function. |
| `.jointDrawSpec()` *(new)* | Build the transformed mean and covariance `Σ_T = B Σ Bᵀ`, and its Cholesky factor, once per run. |
| `.drawJoint()` *(new)* | One joint `(θ, φ)` draw; back-transform each block into a full OMEGA. |
| `.drawableOmegaBlocks()` *(new)* | Apply the coverage policy: return only blocks whose every declared entry is present in `fit$cov` and unfixed, plus the reason each excluded block was held fixed. Shared by both modes. |
| `.omegaWishartSpec()` *(new)* | `"independent_iw"` only. Per block: compute `nu_i` per diagonal, return `nu = min(nu_i)` with validation. |
| `.drawOmega()` *(new)* | `"independent_iw"` only. Delegates the draw to `rxode2::cvPost(nu, omega * (nu - p - 1)/nu)` — the pre-scaling cancels `cvPost`'s `Psi = nu * omega` convention, yielding a mean-centred draw. Scatters results into a full OMEGA template. |
| `.resolveCovarianceParameterSets()` | Dispatch on `control$covarianceDraw`; extend `parameterPartition`. |
| `runSSEControl()` | New `covarianceDraw = c("independent_iw", "joint")` argument, validated by `match.arg()` and rejected when `parameterSource != "covariance"`. New `omegaRseWarn = 0.5` threshold for the weak-identification warning. |
| `runSSE()` / `.validateResumeRequest()` | Record `covarianceDraw` in `run_info`; validate it on resume and `addModels`. A covariance directory with NO recorded value predates OMEGA drawing, so resume ABORTS and requires `restart = TRUE`; `addModels` continues but keeps the recorded value `NA`. |

**Do not hand-roll the inverse-Wishart.** `rxode2::cvPost()` is exported and
already implements it (`type = "invWishart"`, plus `"lkj"` and `"separation"`).
Its convention is `Psi = nu * omega`, so `E[Ω*] = nu/(nu - p - 1) · Ω` — *not*
mean-centred. Verified empirically: at `nu = 20, p = 2` the observed ratio is
1.179 against a predicted 1.176; at `nu = 200`, 1.014 against 1.015. Pre-scaling
the input by `(nu - p - 1)/nu` cancels this exactly, recovering `E[Ω*] = Ω̂` and
the textbook `Var(Ω_ii) = 2Ω̂_ii²/(nu - p - 3)` that the moment-match formula
assumes. Verified: `E[Ω] = 0.2998/0.1200` against a target of `0.30/0.12`, with
all draws positive-definite.

OMEGA name construction (`om.<eta>`, `cov.<eta_i>.<eta_j>`) is generated
locally from `fit$omega` dimnames and intersected against `fit$cov`, rather
than parsing `fit$cov` names directly.

## Error Handling

| Condition | Behaviour |
| --- | --- |
| Weak identification of an OMEGA block | `nu = p + 3 + 2(Ω̂ᵢᵢ/SEᵢ)²` is always `> p + 3` by construction, so a `nu <= p + 3` check catches only underflow or invalid input — not ordinary weak identification. Use an explicit **relative-SE policy** instead: warn when any drawn OMEGA diagonal has relative SE above `omegaRseWarn` (default 0.5), naming the block and eta and reporting the implied behaviour (for `"joint"`, the exact expected mean inflation for a 1×1 block or a scalar-RSE bias proxy for larger blocks; for `"independent_iw"`, the resulting `nu`). Still abort on non-finite `nu`. |
| Fitted OMEGA block not positive-definite | Validate **before** drawing, and abort naming the block and its eta names. Applies to both modes: `"independent_iw"` needs it for `cvPost()`, `"joint"` needs it for `chol()`. |
| A declared OMEGA entry is fixed, or absent from `fit$cov` | **Hold the entire containing block** at its fitted values and record it in `parameterPartition$fixed` with the reason. Decided once by `.drawableOmegaBlocks()`; both modes consume that result and neither re-derives drawability. This supersedes any per-element rule — a single missing off-diagonal holds the whole correlated block, it does not merely drop one eta from a calculation. |
| OMEGA SE zero, `NA`, or non-finite on an otherwise covered block | The block is covered, so it is drawable, but `nu` cannot be moment-matched from that diagonal. Exclude that diagonal from `min(nu_i)`; if no diagonal in the block yields a usable `nu`, keep the block at its fitted values. This is a *distinct* condition from coverage — coverage asks whether the entry exists, this asks whether its SE is numerically usable. |
| `fit$cov` absent or not positive-definite | Existing `.validateCovarianceFit()` behaviour, unchanged. |
| `fit$cov` OMEGA names present that match no eta in `fit$omega` | Abort with the unmatched names, rather than guessing — guards against a future `nlmixr2est` naming change. |
| `"joint"`: `chol()` of an OMEGA block fails (block not positive-definite) | Abort naming the block — a fitted OMEGA that is not PD cannot be log-Cholesky transformed. |
| `"joint"`: transformed covariance `Σ_T` not positive-definite | Abort. Indicates a rank-deficient Jacobian or a degenerate `fit$cov` sub-block. |
| `"joint"`: a block has no `fit$cov` coverage at all | Exclude it from the transform and hold it at its fitted value, recorded in `parameterPartition$fixed`. |
| `covarianceDraw` supplied when `parameterSource != "covariance"` | Abort from `runSSEControl()`, matching the existing gating of `inFilter`/`randomEstimationInits` to `"rawres"`. |

## Testing

Unit tests, using constructed fits (no live estimation):

- `.alignedCovariance()` partitions a mixed theta/OMEGA `fit$cov` correctly,
  and no longer aborts on `om.`/`cov.` names.
- Block detection recovers 1×1 blocks, a single 2×2 block, and mixed cases.
- `nu` matches the closed form for known inputs; `min` selection picks the
  binding element.
- Fixed etas are excluded and restored at fitted values.
- Each documented error condition raises with the offending name.

- `"joint"`: the log-Cholesky transform round-trips (`.phiToOmega(.omegaToPhi(Ω)) == Ω`)
  for a 1×1 block, a 2×2 block, and a 3×3 block with off-diagonals.
- `"joint"`: the numeric Jacobian matches a finite-difference check at a
  coarser step, confirming the step size is sane.
- `covarianceDraw` defaults to `"independent_iw"`, accepts `"joint"`, rejects
  anything else, and is refused when `parameterSource != "covariance"`.

Required edge cases — each has a specific failure mode this design must not
regress into:

- **Theta-only `fit$cov`.** Thetas must still be drawn and all OMEGA recorded
  fixed. The existing `fake_sse_fit()` fixture already has a theta-only `cov`,
  so the current test at `test-run-sse.R:168` covers this and must keep
  passing.
- **OMEGA-only covariance, zero drawn thetas.** Guards the
  `drawn[-seq_len(nTheta)]` bug: with `nTheta == 0`, `-seq_len(0)` is
  `-integer(0)` and R returns an **empty** vector, silently discarding the
  entire OMEGA draw.
- **A correlated block containing a fixed eta variance.** Must hold the whole
  block fixed, not abort and not partially draw.
- **A declared block with entries missing from `fit$cov`.** Same policy as
  above.
- **1×1, 2×2 and 3×3 blocks**, confirming element ordering.
- **Near-boundary and ill-conditioned OMEGA**, and OMEGA scales spanning
  roughly `1e-8` to `1e2`, exercising the scale-aware Jacobian step.
- **Accurate drawn/fixed partition reporting** in every case above.
- **Resume and `addModels` with a changed `covarianceDraw`.**

Statistical tests (seeded, tolerance-based):

- `"independent_iw"`: over many draws from a known `Omega0` and target SEs, empirical
  `E[Omega]` recovers `Omega0` and the binding element's empirical SD recovers
  its target SE.
- `"joint"`: over many draws from a covariance with a **deliberately imposed
  THETA↔OMEGA correlation**, the recovered correlation matches the imposed one
  within tolerance. This is the property the whole mode exists for, so it must
  be asserted directly rather than inferred.
- `"joint"`: **every** draw is positive-definite. Assert across a few thousand
  draws — this is the guarantee that removes rejection sampling, so it should
  fail loudly if the transform is ever broken.

Both use a fixed seed and a tolerance wide enough to be robust across
platforms.

Integration test:

- `runSSE(parameterSource = "covariance")` completes end to end on a real fit
  under **both** modes, and the recorded `parameterPartition` lists the
  expected drawn OMEGA entries.

### Existing assertions and documentation that become false

Once OMEGA is drawn, several current statements stop being true. Precisely:

| Location | Status |
| --- | --- |
| `tests/testthat/test-run-sse.R:888` — asserts drawn OMEGA has exactly one unique value | Now false, but **currently masked**: it sits inside the end-to-end test that `skip()`s because of the `nlmixr2utils` `cli` defect. It must still be corrected, or it becomes a trap the moment that defect is fixed. |
| `tests/testthat/test-run-sse.R:168` — "draws theta and keeps omega fixed" | **Stays correct.** Its `fake_sse_fit()` fixture has a theta-only `cov`, so the coverage policy keeps OMEGA fixed. Retain it as the theta-only regression test. |
| `R/run-sse.R:291` — cli message "OMEGA and SIGMA stay at the fitted point estimates" | Must be rewritten to describe the selected mode. |
| `README.md:97` — "OMEGA and SIGMA stay fixed" | Must be rewritten. |
| `vignettes/runSSE.Rmd:203, 219` — "theta-only uncertainty" and "For full uncertainty in OMEGA and SIGMA as well, use …" | Must be rewritten. |

The plan must therefore **not** claim the suite retains an unchanged
expected-failure list. The currently-executing set happens to be unaffected —
only because the false assertion is masked by a skip — and the documentation
claims are false regardless.

## Risks

- **Naming-convention dependency.** `om.<eta>` / `cov.<eta_i>.<eta_j>` is
  produced by unexported `@noRd` helpers in `nlmixr2est`
  (`.foceiOmegaCovNames()`), not a documented public contract. Mitigated by
  generating expected names locally and erroring clearly on any unmatched name
  rather than misparsing silently.
- **`"joint"` inflates OMEGA means, and not slightly.** The log-Cholesky
  back-transform is non-linear, so `E[Ω*] ≠ Ω̂`. For a 1×1 block, `φ = ½ log Ω`
  and the centred Gaussian draw makes `Ω*` lognormal with log-sd `SE/Ω̂`, so:

  ```
  E[Ω*] / Ω̂  =  exp( SE² / (2 Ω̂²) )
  ```

  | OMEGA relative SE | mean inflation |
  | ---: | ---: |
  | 20% | 2.0% |
  | 50% | 13.3% |
  | 100% | 64.9% |

  An earlier draft quoted "about +2%" as if it were general; that was the
  20%-RSE case only. Poorly identified variance components routinely reach
  50–100% relative SE, and those are precisely the fits where propagating
  parameter uncertainty matters most. This is why `"joint"` is opt-in and must
  warn when relative SE is large.
- **Delta-method approximation.** `Σ_T` is a first-order local approximation.
  After the non-linear back-transform, `"joint"` does not exactly reproduce
  raw-scale means, variances, or Pearson correlations, and does not reproduce
  the fitted asymptotic distribution. It *incorporates* the cross-covariance;
  it does not *preserve* it. The approximation degrades as OMEGA uncertainty
  grows — the same regime as the mean inflation above.
- **`"independent_iw"` is mean-centred, not "unbiased".** `E[Ω*] = Ω̂` is a
  property of the generated distribution, not an estimator property over
  repeated datasets.
- **Single `nu` per block** (`"independent_iw"` only) ignores off-diagonal
  OMEGA standard errors and the covariance among OMEGA estimators, and lets
  the least precise diagonal inflate the whole block. `min(nu_i)` ensures only
  that each constructed diagonal marginal SD *matches or exceeds* its
  requested value — it is not guaranteed conservative for non-linear
  predictions, correlations, ΔOFV, Type-I error, or power.
- **Reproducibility break.** Same-seed runs will not reproduce prior releases'
  draws. Requires a NEWS entry. Note the two modes also differ from each other
  at the same seed, which is expected.
- **Poorly identified variance components** now surface as warnings where they
  were previously ignored. Intended, but a visible behaviour change.
- **THETA draws remain unconstrained, in both modes.** The existing covariance
  mode already draws thetas from an unconstrained multivariate Normal, and
  this design does not change that. Bounded thetas and residual-error standard
  deviations can therefore be drawn outside their valid domains — a negative
  `add.sd`, or a theta below its declared `lower`. Extending the mode to "full
  parameter uncertainty" makes this pre-existing limitation more prominent, so
  it must be stated in the user documentation. Transformed or truncated theta
  draws are a later design.
- **`covarianceDraw` is part of the reproducibility contract.** Switching modes
  changes every simulated parameter set at the same seed, so the mode must be
  recorded in `run_info` and validated on resume and `addModels` — the same
  treatment `parameterSource`, `randomEstimationInits`, and `updateFix` already
  receive, and the same treatment `rxThreads` received for the same reason.
  Runs completed before `covarianceDraw` existed carry no recorded value, and
  a covariance run with no recorded value **necessarily held OMEGA fixed** —
  it predates OMEGA drawing entirely. Resuming one would append OMEGA-varying
  replicates to OMEGA-fixed replicates, putting two simulation distributions
  in a single study. So the missing value is not merely "unknown": **abort the
  resume and require `restart = TRUE`.**

  This differs deliberately from `rxThreads`, where a missing record warns
  rather than aborts. There the thread count is unknown but the *kind* of draw
  is unchanged; here the previous behaviour is known, and known to be
  incompatible.

  `addModels` is safe to continue, because it refits against saved datasets and
  simulates nothing new — but the recorded `covarianceDraw` must stay `NA` for
  such a run rather than being rewritten to the current mode, which would
  falsely describe how those datasets were produced.

## Out of Scope

- The `nlmixr2utils:::.abortRawResults()` `cli` `.envir` defect masking the
  covariance/rawres integration test. File separately.
- `runSSEControl(initialEtas = ...)`, still hard-disabled — the analogue of
  PsN's `-initial_etas`. Unrelated to parameter uncertainty; separate work.
- Repairing the legacy `fit$sigma` accessor in `nlmixr2est`.
- A bias correction for `"joint"`'s non-linear back-transform. Quantified and
  documented as a known limitation instead; revisit only if it proves to
  matter in practice.
- Alternative unconstrained parameterizations (matrix-log, separation
  strategy). Log-Cholesky is sufficient and is the cheapest to compute.
