# SSE Initial ETAs Design

## Goal

Enable `runSSEControl(initialEtas = TRUE)` so the per-subject ETA values used to
simulate each replicate are handed to the estimation step as starting values.
This is the analogue of PsN's `sse -initial_etas`, and the control argument
already exists in `nlmixr2sse` but is hard-disabled.

This is a companion to `2026-08-28-sse-parameter-uncertainty-design.md`. The two
are independent: this feature concerns the simulation-to-estimation handoff and
applies to every `parameterSource` mode.

## Background

`runSSEControl()` accepts `initialEtas` and then refuses it outright:

```r
if (isTRUE(initialEtas)) {
  .abortSSE(
    "{.arg initialEtas} is still reserved and must remain FALSE in the current release."
  )
}
```

The documentation describes it as "reserved ... pending stability work". The
supporting pieces are now available on both sides.

### What PsN does

`sse -initial_etas` tables the simulated per-subject ETAs to a `.simphi` file
and attaches them to every estimation model with an `$ETAS` record
(`lib/tool/sse.pm:1024`, `:1285`):

```perl
$orig_est_models[$j]->problems->[0]->set_records(
    type => 'etas',
    record_strings => [ "FILE=mc-" . ($j + 1) . ".simphi" ] );
```

PsN attaches this unconditionally to the simulation model and to every
alternative, without checking that their random-effect structures match.

### What nlmixr2 provides

`foceiControl(etaMat = )` takes an initial ETA matrix — one row per subject,
one column per eta — and `checkmate::assertMatrix(etaMat, mode = "double")`
validates it. `fit$etaMat` has exactly this shape.

The simulated ETAs are recoverable from the `rxSolve()` result: the solved
object carries a `params` data frame with one row per subject, an `id` column,
and one column per eta:

```
  id       tka      tcl       tv    add.sd     eta.ka     eta.cl
1  1 0.4024209 1.027094 3.429623 0.7819965  0.2802193 -0.1602785
2  2 0.4024209 1.027094 3.429623 0.7819965  0.2408827  0.2280125
```

## Chosen Approach

Capture the simulated ETAs during `.simulationRecord()`, carry them on the
simulation record, and pass them into the estimation step in
`.fitTaskRecord()` via `foceiControl(etaMat = )`.

### Capturing the ETAs

`.simulationRecord()` currently discards the information. Its `rxSolve()` call
requests `returnType = "data.frame"` and immediately wraps the result in
`as.data.frame()`, which keeps only the solved rows — the `params` component is
lost. The call must retain the `rxSolve` object long enough to read
`solved$params`, then convert as before.

### Subject alignment

`solved$params$id` is a **factor whose labels are the original subject IDs**,
in sorted order — not a positional index. Verified with non-sequential IDs
(`101, 205, ..., 1300`): the factor labels carry those values.

Rows must therefore be aligned by `as.character(solved$params$id)` against the
estimation data's own ID column. Using `as.integer()` on the factor would
silently yield level indices (`1..n`) and misalign every subject whenever IDs
are not `1..n` — a correctness failure that would not error, only bias results.
The implementation must convert via character, and a test must cover
non-sequential IDs specifically.

### Column alignment

`etaMat` columns must match, in order, the etas of the model being fitted.
Model etas are the `fit$iniDf` rows where `neta1 == neta2`.

- Model etas ⊆ simulated etas → subset and reorder the captured matrix to the
  model's own eta order, then pass it.
- Model has an eta absent from the simulation → fill that column with `0`,
  the eta's prior mean and the value it would otherwise start at. The model
  still receives true starting values for every eta the simulation shares with
  it, and is no worse off than the default for the rest.

Filling with zero keeps the feature available for a legitimate SSE use case:
testing the power to detect a *new* random effect, where the alternative model
deliberately adds an eta the simulation model does not have. Aborting instead
would make `initialEtas` unusable for exactly that comparison.

Zero-filling is reported, not silent: the run records which etas were
zero-filled for which model (see Reporting below), so a reader can tell that
such a model started from a partly-true, partly-default point.

### Estimation-method gating

`etaMat` is a `foceiControl` argument. `saemControl()` has no equivalent, and
neither do the other non-FOCEI estimators.

If `initialEtas = TRUE` and any scheduled model uses a method without `etaMat`
support, abort **before** any simulation work begins, naming the model and its
method. Silently starting some models from the true ETAs and others from zero
would give them different starting conditions while their OFVs are compared
against each other, quietly biasing the comparison the tool exists to make.

The check runs as a pre-flight pass over every scheduled model — the simulation
model and every alternative — in `runSSE()`, before the first replicate is
simulated and before any run directory work. Requirements for the message:

- It names **every** unsupported model and its method in one error, not just
  the first encountered, so a user with several alternatives fixes them in one
  pass rather than one abort at a time.
- It names the supported methods explicitly.
- It states the two ways forward: change the model's estimation method, or set
  `initialEtas = FALSE`.

The supported-method set is derived by checking for an `etaMat` formal on the
method's control constructor, rather than hard-coding a list of method names,
so a method that gains `etaMat` support in a future `nlmixr2est` is picked up
without a change here. If that introspection cannot resolve a method's control
constructor, the method is treated as unsupported — failing closed, since the
cost of a wrong "supported" answer is a silently biased run.

### Avoiding a spurious warning

`foceiControl()` warns when `etaMat` is supplied without an explicit
`maxInnerIterations` (`foceiControl.R:1750`). The implementation passes
`maxInnerIterations` explicitly, preserving any value the user already set on
the model's control and otherwise using the `foceiControl` default, so the
warning does not fire on every replicate.

## Component Design

| Function | Change |
| --- | --- |
| `runSSEControl()` | Remove the hard block on `initialEtas`. Validate it as a flag. |
| `.simulationRecord()` | Retain the `rxSolve` object, read `solved$params`, extract eta columns aligned by `as.character(id)`, store as `simRecord$etaMat`. |
| `.etaMatForModel()` *(new)* | Subset and reorder a captured eta matrix to one model's eta names; zero-fill etas the simulation lacks and report which were filled. |
| `.assertInitialEtasSupported()` *(new)* | Pre-flight check that every scheduled model's estimation method supports `etaMat`. |
| `.fitTaskRecord()` | When `initialEtas` is set, inject `etaMat` and an explicit `maxInnerIterations` into the model's control. |
| `runSSE()` | Run the pre-flight check before the simulation loop. |

`initialEtas` is orthogonal to `parameterSource` and composes with `fixed`,
`rawres`, and `covariance` alike; no interaction handling is required.

### Reporting

`runInfo` records that `initialEtas` was active, and for each model the etas
that were supplied from the simulation versus zero-filled. Because a
zero-filled model starts from a partly-true point while a fully-matched model
starts from an entirely true one, and their OFVs are then compared, this
difference must be visible in the run record rather than inferred.

### Resume and add-models interaction

`simRecord$etaMat` becomes part of the cached simulation record, so resumed
runs reuse the captured ETAs. `addModels = TRUE` fits new alternatives against
saved datasets: if those records predate this feature they carry no `etaMat`,
so an add-models run requesting `initialEtas = TRUE` against such a directory
must abort with a clear message rather than silently starting from zero.

## Error Handling

| Condition | Behaviour |
| --- | --- |
| Any scheduled model uses an estimation method without `etaMat` support | Abort before simulating, naming the model and method. |
| A model has etas absent from the simulation model | Zero-fill those columns and record them; not an error. |
| `solved$params` missing, or lacks an `id` column | Abort — the rxode2 result shape is not as expected. |
| Simulated subject IDs do not match the estimation data's IDs | Abort reporting the mismatch, rather than aligning positionally. |
| `addModels = TRUE` against saved records lacking `etaMat` | Abort explaining that the original run did not capture ETAs. |
| Model has no etas at all | Skip `etaMat` for that model; not an error. |

## Testing

Unit tests, with constructed inputs:

- Eta capture pulls the correct columns from a `params`-shaped data frame.
- **Non-sequential IDs** (e.g. `101, 205, ..., 1300`) align correctly; a test
  asserts that positional/`as.integer()` alignment would have produced a
  different, wrong matrix.
- Column subsetting reorders to a model's eta order, and drops etas the model
  does not have.
- A model with an eta absent from the simulation gets a zero-filled column in
  the right position, with the shared etas still carrying their true values,
  and the zero-fill is recorded.
- Each documented abort fires with the offending name.
- The method pre-flight names **all** unsupported models in one error, not
  only the first, and fires before any simulation work.
- A method whose control constructor cannot be resolved is treated as
  unsupported.
- A model with no etas is skipped without error.

Integration tests:

- `runSSE(..., control = runSSEControl(initialEtas = TRUE))` completes end to
  end on a real FOCEI fit with a matching alternative.
- The same call with a `saem` model aborts before simulating.
- `initialEtas = TRUE` composes with `parameterSource = "covariance"`.

## Risks

- **Optimizer starting-point bias.** Starting the estimation at the true
  simulated ETAs is a deliberate departure from a naive analysis. It can
  improve convergence, but it can also flatter a model relative to real-world
  use where true ETAs are unknown. This is inherent to the feature and to PsN's
  equivalent; it belongs in the user documentation, not in the code.
- **Unspecified reason for the original hard-disable.** The control
  documentation cites "stability work" without detail. The design may
  therefore be re-encountering a known problem. Mitigated by the pre-flight
  gating and the alignment tests, but worth watching during implementation; if
  a concrete instability surfaces, revisit rather than force it through.
- **`solved$params` is not a documented public contract.** It is a component of
  the `rxSolve` result rather than a documented accessor. Mitigated by
  asserting its shape and aborting clearly if absent.

## Out of Scope

- Parameter uncertainty in THETA/OMEGA — see
  `2026-08-28-sse-parameter-uncertainty-design.md`.
- Adding `etaMat` support to `saemControl()` in `nlmixr2est`.
- PsN's `-keep_tables` / `special_table` options, unrelated to ETAs.
