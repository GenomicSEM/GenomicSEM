# lavaan_rust experiment plan

## Goal

Create an experimental GenomicSEM branch that measures how far we can get by
replacing only the lavaan layer first.

The initial rule is deliberately strict:

- GenomicSEM functions keep their current logic.
- New wrappers such as `commonfactor_rust()` and `usermodel_rust()` should be
  line-for-line equivalents of the existing functions except that they call
  `lavaan_rust` entrypoints such as `sem_rust()` and `lavInspect_rust()`.
- `_rust()` wrappers are strict: if execution reaches a lavaan feature that the
  Rust backend does not yet implement, they should error instead of silently
  delegating back to lavaan.
- `lavaan_rust` should implement only the subset of lavaan that GenomicSEM
  actually uses, with an R-facing API that mirrors lavaan closely enough for
  those wrappers to stay mechanically comparable.
- Broader GenomicSEM rewrites, batching, and algorithmic changes are deferred
  until we have measured the isolated benefit of a faster lavaan-compatible
  backend.

This gives us a clean experiment: if `f_rust()` is faster than `f()` while the
wrapper logic is otherwise unchanged, the measured improvement belongs to the
lavaan replacement rather than to changed GenomicSEM orchestration.

## Current lavaan surface in upstream GenomicSEM

The current upstream code imports all of `lavaan`, but the observed usage surface
is much smaller:

| API | Static call count | Notes |
| --- | ---: | --- |
| `lavInspect()` | 53 | Most important accessor surface |
| `sem()` | 50 | Main fitting entrypoint |
| `inspect()` | 29 | Mostly parameter-list and SE extraction |
| `parTable()` | 15 | Used both for reporting and model mutation |
| `lavaan()` | 4 | Refit path using a precompiled base model |
| `fitted()` | 4 | Implied covariance matrices |
| `lav_model_get_parameters()` | 5 | Defined-parameter SE support |
| `lav_func_jacobian_complex()` | 5 | Defined-parameter SE support |
| `resid()` | 1 | Common-factor GWAS Q model |

The call sites are concentrated in:

| File | Static call count |
| --- | ---: |
| `R/rgmodel.R` | 47 |
| `R/usermodel.R` | 34 |
| `R/commonfactor.R` | 27 |
| `R/enrich.R` | 18 |
| `R/commonfactorGWAS_main.R` | 18 |
| `R/userGWAS_main.R` | 11 |
| `R/userGWAS.R` | 9 |

### Observed accessor keys and object expectations

`lavInspect()` currently needs to support at least:

- `"cor.lv"`
- `"delta"`
- `"WLS.V"`
- `"converged"`
- `"fit"`
- `"vcov"`

`inspect()` currently needs to support at least:

- `"list"`
- `"se"`
- default `inspect(fit)[[1]]` covariance ordering behavior

The wrappers also depend on behavior beyond ordinary function calls:

- `class(fit)[1] == "lavaan"`
- direct access to slots such as `fit@Model`, `fit@Options`, `fit@ParTable`,
  `fit@Data`, and `fit@Model@def.function`
- reusing a precompiled base model through `lavaan(..., slotOptions=...,
  slotParTable=..., slotData=..., slotModel=...)`
- mutating `parTable()` output and passing the mutated table back into `sem()`
- user-defined parameters with `:=`
- inequality constraints such as residual-variance lower bounds

That means this is not just an optimizer project. The compatibility layer is part
of the backend contract.

## Proposed repository shape

Start with a nested R package so the experiment is isolated and easy to remove.
The conceptual backend is `lavaan_rust`, but the R package is named
`lavaanrust` because R package names cannot contain underscores:

```text
lavaanrust/
  DESCRIPTION
  NAMESPACE
  R/
    wrappers.R
    classes.R
  src/rust/
    Cargo.toml
    src/
      lib.rs
      parser.rs
      model.rs
      fit.rs
      inspect.rs
      par_table.rs
```

The R package should use `extendr` and expose R-callable functions with lavaan-like
signatures but `_rust` suffixes:

- `sem_rust()`
- `lavaan_rust()`
- `lavInspect_rust()`
- `inspect_rust()`
- `parTable_rust()`
- `fitted_rust()`
- `resid_rust()`
- `lav_model_get_parameters_rust()`
- `lav_func_jacobian_complex_rust()`

The result object should initially be an S4 compatibility class with the minimum
slots that GenomicSEM observes. It does not need to be a general-purpose lavaan
replacement; it needs to behave like lavaan for this package's supported slice.

## Staged roadmap

### Phase 0: freeze the compatibility contract

1. Add golden tests around the current lavaan-backed outputs before changing
   behavior.
2. Record representative fixtures for:
   - simple common-factor models
   - `usermodel()` with and without `std.lv`
   - defined parameters via `:=`
   - residual lower-bound retry models
   - common-factor GWAS Q-model refits
   - base-model reuse through `lavaan(...)`
3. Snapshot every accessor surface GenomicSEM consumes:
   - parameter tables
   - inspect-list rows and column ordering
   - delta matrices
   - WLS weight matrices
   - fit summaries
   - implied covariance matrices
   - convergence flags

### Phase 1: scaffold `lavaan_rust`

1. Add the nested R package and `extendr` build plumbing.
2. Define the compatibility S4 fit object and R wrappers.
3. Add the first narrow tests that prove R can construct, inspect, serialize,
   and compare rust-backed fit objects.

### Phase 2: implement the first useful SEM slice

Implement enough syntax and fitting support for covariance-only DWLS models used
by `commonfactor()`:

- operators: `=~`, `~~`
- fixed coefficients such as `1*F1`
- free loadings and residual variances
- latent variance normalization
- covariance-only input via `sample.cov`
- `WLS.V`
- `std.lv`
- convergence metadata

Required output support for this slice:

- `sem_rust()`
- `lavInspect_rust(..., "delta")`
- `lavInspect_rust(..., "WLS.V")`
- `lavInspect_rust(..., "converged")`
- `lavInspect_rust(..., "cor.lv")`
- `parTable_rust()`
- `inspect_rust(..., "list")`
- `fitted_rust()`

The first GenomicSEM wrapper should be `commonfactor_rust()`. It is the cleanest
place to test the idea because it exercises real GenomicSEM logic while avoiding
the widest part of user-defined model syntax.

### Phase 3: make model reuse and mutation work

Add the subset needed by GWAS loops:

- `lavaan_rust()` refits from a precompiled model object
- stable `ParTable` reuse
- mutating a parameter table and passing it back into `sem_rust()`
- `resid_rust()`

This unlocks `commonfactorGWAS_rust()` and gives the experiment a realistic
high-throughput workload without changing loop structure yet.

### Phase 4: cover user-defined models

Broaden syntax and compatibility for `usermodel_rust()` and `userGWAS_rust()`:

- regressions via `~`
- labels and equality constraints
- inequality constraints
- defined parameters via `:=` for the supported arithmetic subset
- `lav_model_get_parameters_rust()`
- `lav_func_jacobian_complex_rust()`
- `inspect_rust(..., "se")`
- `lavInspect_rust(..., "fit")`
- `lavInspect_rust(..., "vcov")`

This is the first phase where parser breadth and lavaan-internal compatibility
become the dominant risk.

The compiler-first refinement for this phase is tracked separately in
`LAVAAN_FAST_PLAN.md`: before adding more one-off model families, lower
lavaan-style parameter tables into a reusable intermediate representation and
move existing supported fits onto that representation incrementally.

### Phase 5: broaden only after measurements

Once the isolated swap experiment is measured:

- Decide whether the remaining GenomicSEM functions are worth cloning into
  `_rust` variants.
- Compare isolated backend gains against gains from batching and orchestration
  changes.
- Only then decide whether to replace original functions, keep an alternate
  backend, or continue with deeper GenomicSEM rewrites.

## Test strategy

### Compatibility tests

For every supported wrapper pair:

```r
old <- f(...)
new <- f_rust(...)
```

Compare:

- class and names
- parameter-row ordering
- warning/error behavior where practical
- exact string fields
- matrix dimensions and dimnames
- numeric values

Use byte-identical snapshots for textual and structural outputs where possible.
Use strict numeric tolerances only where floating-point implementation details
make byte equality unrealistic, and record those tolerances explicitly.

### Differential tests

- Generate small valid lavaan model strings within the supported syntax subset.
- Fit each model through lavaan and `lavaan_rust`.
- Compare parsed parameter tables, implied covariance matrices, gradients, and
  final estimates.
- Add negative tests for unsupported syntax so the rust backend fails loudly
  instead of silently accepting a different model.

### Fuzzing

- Fuzz parser input for supported operators, labels, comments, whitespace, and
  invalid token sequences.
- Fuzz small positive-definite covariance matrices and WLS diagonal inputs.
- Differentially fuzz small identified models against lavaan whenever both
  backends should accept the same syntax.

### Performance tracking

Keep a timestamped benchmark log from the start. The initial comparison set
should contain:

- lavaan-backed `f()`
- isolated rust-backed `f_rust()`
- later, if desired, any higher-level optimized GenomicSEM path

For this branch, keep the first benchmarks single-fit and same-control-flow. The
point is to isolate backend speed before layering in batching or multithreading.

## First implementation packet

The smallest packet that is both technically useful and experimentally clean is:

1. Add the nested `lavaan_rust` package skeleton with `extendr`.
2. Implement a compatibility S4 fit object plus empty accessor wrappers.
3. Add golden fixtures from upstream lavaan for simple one-factor DWLS models.
4. Implement `sem_rust()` for the narrow `commonfactor()` slice.
5. Add `commonfactor_rust()` as an otherwise-identical wrapper around
   `commonfactor()`.
6. Benchmark `commonfactor()` vs `commonfactor_rust()` before adding any
   GenomicSEM-specific batching or orchestration changes.

That packet answers the first real question cleanly: how much speed can we buy
from a lavaan-compatible rust backend alone?
