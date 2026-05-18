# lavaan_rust progress log

## 2026-05-04 12:34 PDT

### Implemented

- Added the nested experimental R package `lavaanrust/`.
- Added an `extendr` Rust backend with `fit_one_factor_dwls()`.
- Added a minimal S4 compatibility object, `lavaan_rust_fit`, and the first
  suffixed compatibility surface:
  - `sem_rust()`
  - `lavaan_rust()`
  - `lavInspect_rust()`
  - `inspect_rust()`
  - `parTable_rust()`
  - `fitted_rust()`
  - `resid_rust()`
  - defined-parameter helper stubs that currently strict-error
- Added `commonfactor_rust()` without copying the original function body:
  it clones the existing `commonfactor()` closure and rebinds only the lavaan
  symbols in a child environment.
- Added a second native closed-form DWLS slice for observed-covariance models,
  covering the common-factor CFI/null model and its parameter-table refit.
- Switched `_rust()` policy to strict mode: unsupported paths now error instead
  of silently delegating to upstream lavaan.
- Added frozen one-factor fixtures from lavaan 0.6-21 and test coverage for:
  - the rust-backed one-factor DWLS slice
  - strict errors for unsupported syntax

### Current support boundary

`sem_rust()` currently owns:

- one-factor covariance-only DWLS models used by the primary and standardized
  `commonfactor()` fits
- the observed-covariance DWLS family used by the common-factor CFI/null model
  and its parameter-table refit

Unsupported syntax now errors deliberately; `_rust()` wrappers are meant to
measure the Rust surface that actually exists, not hide missing coverage by
delegating back to lavaan.

### Validation

- `lavaanrust/tests/testthat`: 12 passing tests.
- `R CMD build lavaanrust` followed by
  `R CMD check lavaanrust_0.0.0.9000.tar.gz --no-manual`: `Status: OK`.
- Synthetic `commonfactor()` smoke comparison:
  - max unstandardized estimate difference: `3.74e-09`
  - max sandwich SE difference: `3.72e-09`

### Local benchmarks

Synthetic 3-trait fixture, macOS laptop, 2026-05-04:

| Benchmark | lavaan | rust-backed | Speedup |
| --- | ---: | ---: | ---: |
| direct one-factor `sem()` fit | `0.031630 s` | `0.000828 s` | `38.20x` |
| end-to-end `commonfactor()` smoke, first packet | `0.164460 s` | `0.064800 s` | `2.54x` |
| direct one-factor `sem()` fit after null-model slice | `0.017726 s` | `0.000522 s` | `33.96x` |
| end-to-end `commonfactor()` smoke after null-model slice | `0.092230 s` | `0.006320 s` | `14.59x` |

The direct-fit gain is the clean measurement of the Rust one-factor backend.
The second end-to-end benchmark shows the impact of removing the remaining
lavaan work from `commonfactor_rust()` itself.

### Next packet

1. Add strict base-model reuse support through `lavaan_rust()` for GWAS loops.
2. Broaden `sem_rust()` from current covariance-only families toward the
   `commonfactorGWAS()` model slice.
3. Add `commonfactorGWAS_rust()` once the model-reuse path is native.

## 2026-05-04 12:59 PDT

### Implemented

- Added a strict native `commonfactorGWAS()` backend family:
  - marker-scaled factor/SNP DWLS first-stage fit
  - slot-style base-model reuse through `lavaan_rust()`
  - Q-model refit with direct SNP effects and trait residual variances free
- Extended `lavaan_rust_fit` with the compatibility slots needed by the current
  GenomicSEM reuse path: `Options`, `Data`, and `Model`.
- Added `commonfactorGWAS_rust()` using the original GenomicSEM implementation
  with rust-bound helper closures for `.commonfactorGWAS_main()` and
  `.rearrange()`.
- Kept the wrapper boundary strict: `commonfactorGWAS_rust()` currently requires
  `parallel = FALSE`; unsupported parallel execution errors instead of drifting
  into the original worker path.
- Added a narrow `class()` compatibility shim inside the rust-bound wrapper
  environment so the unchanged upstream `class(fit)[1] == "lavaan"` checks see
  rust fit objects as successful lavaan-compatible fits.

### Current support boundary

`sem_rust()` now owns:

- one-factor covariance-only DWLS models used by `commonfactor()`
- the observed-covariance DWLS family used by the common-factor CFI/null model
  and its parameter-table refit
- the marker-scaled one-factor SNP model used by `commonfactorGWAS()`
- the `commonfactorGWAS()` Q-model parameter-table refit

`lavaan_rust()` now owns the strict slot-reuse path for the supported
`commonfactorGWAS()` DWLS base model. Unsupported reuse patterns still error and
do not delegate to lavaan.

`commonfactorGWAS_rust()` currently supports the sequential path only. The
existing `foreach` worker setup needs an explicit rust-bound export strategy
before the parallel path can be considered supported.

### Validation

- `lavaanrust/tests/testthat`: 20 passing tests.
- Frozen lavaan 0.6-21 `commonfactorGWAS()` fixture:
  - first-stage parameter estimates match within `2e-06`
  - Q-model refit parameter estimates match within `2e-06`
- Synthetic sequential GenomicSEM smoke:
  - `commonfactorGWAS()` vs `commonfactorGWAS_rust()` factor-SNP estimate max
    absolute difference: `1.57e-08`
  - Q-statistic absolute difference: `4.82e-08`

### Local benchmarks

Synthetic 3-trait fixture, macOS laptop, 2026-05-04:

| Benchmark | lavaan | rust-backed | Speedup |
| --- | ---: | ---: | ---: |
| direct `commonfactorGWAS()` first-stage `sem()` fit | `0.018000 s` | `0.001000 s` | `18.00x` |
| end-to-end sequential `commonfactorGWAS()` smoke, 1 SNP | `0.067000 s` | `0.005000 s` | `13.40x` |

### Next packet

1. Decide whether to push the experiment deeper into `commonfactorGWAS()`
   parallel worker support or move next to `usermodel_rust()`.
2. Broaden parser and constraint support only where another GenomicSEM wrapper
   actually needs it.

## 2026-05-04 13:20 PDT

### Implemented

- Added the first strict `usermodel_rust()` slice while keeping the original
  GenomicSEM `usermodel()` body intact.
- Added support for simple one-factor user models such as
  `F1 =~ A + B + C`:
  - marker scaling when `std.lv = FALSE`
  - latent-variance scaling when `std.lv = TRUE`
- Added the minimal `standardizedSolution_rust()` compatibility surface needed
  by the unchanged `usermodel()` standardized-output path.
- Split the implicit `std.lv = TRUE` constructor from the explicit
  `F1 ~~ 1*F1` constructor so parameter-table row order matches lavaan for both
  syntactic forms.
- Tightened the class shim used by rust-bound wrappers so exact upstream checks
  such as `class(x) != "lavaan"` continue to behave like the original code.

### Current support boundary

`usermodel_rust()` currently supports simple one-factor DWLS models with no:

- labels or equality constraints
- inequality constraints
- user-defined parameters via `:=`
- regressions
- multiple latent factors

Unsupported user-model syntax still errors through `sem_rust()` and does not
fall back to lavaan.

### Validation

- `lavaanrust/tests/testthat`: 28 passing tests.
- Frozen lavaan 0.6-21 fixtures for:
  - marker-scaled one-factor user models
  - `std.lv = TRUE` one-factor user models
  - standardized-solution rows consumed by `usermodel()`
- Synthetic `usermodel()` smoke, simple marker-scaled one-factor model:
  - max unstandardized estimate difference: `7.17e-09`
  - max sandwich SE difference: `1.60e-09`
  - max `STD_Genotype` difference: `1.11e-08`
  - max `STD_All` difference: `5.97e-09`

### Local benchmarks

Synthetic 3-trait fixture, macOS laptop, 2026-05-04:

| Benchmark | lavaan | rust-backed | Speedup |
| --- | ---: | ---: | ---: |
| direct user-model `sem()` fit | `0.017000 s` | `0.001000 s` | `17.00x` |
| end-to-end simple `usermodel()` smoke | `0.067000 s` | `0.005000 s` | `13.40x` |

### Next packet

1. Expand the user-model parser one feature at a time, starting with the
   smallest syntax family that unlocks a real GenomicSEM workflow beyond the
   current one-factor case.
2. Likely next targets are either constrained one-factor models or the
   user-defined-parameter machinery needed for `:=`.

## 2026-05-04 13:30 PDT

### Implemented

- Added the first strict `userGWAS_rust()` slice while preserving the original
  GenomicSEM `userGWAS()` body.
- Reused the existing marker-scaled one-factor SNP optimizer for the
  unrestricted first-stage `F1 =~ ...` plus `F1 ~ SNP` fit.
- Added a new native fixed-measurement DWLS solver for the default
  `fix_measurement = TRUE` path:
  - fixed loadings from the no-SNP measurement model
  - free trait residual variances
  - free latent residual variance
  - free factor-SNP regression
  - free SNP variance
- Added strict native slot reuse for the fixed-measurement base model through
  `lavaan_rust()`, which is the path exercised inside the per-SNP loop.
- Kept the wrapper boundary explicit. `userGWAS_rust()` currently requires:
  - `parallel = FALSE`
  - `fix_measurement = TRUE`
  - `Q_SNP = FALSE`
  - `estimation = "DWLS"`
  - `TWAS = FALSE`

### Current support boundary

The new `userGWAS_rust()` slice supports simple one-factor models such as:

```r
F1 =~ A + B + C
F1 ~ SNP
```

Unsupported user-GWAS paths still error and do not fall back to lavaan. The
current packet intentionally does not yet cover unconstrained measurement
refits, `Q_SNP`, TWAS, parallel workers, multiple latent factors, or more
general user-model syntax.

### Validation

- `lavaanrust/tests/testthat`: 37 passing tests.
- Frozen lavaan 0.6-21 fixtures for:
  - the unrestricted one-factor user-GWAS first stage
  - the fixed-measurement parameter-table refit
  - fixed-measurement slot reuse
- Synthetic sequential GenomicSEM smoke, 1 SNP:
  - max estimate difference: `5.20e-09`
  - max sandwich SE difference: `3.34e-10`
  - max chi-square difference: `2.06e-09`
- `R CMD build lavaanrust` followed by
  `R CMD check lavaanrust_0.0.0.9000.tar.gz --no-manual`: clean after adding
  the generated Rd page for the new exported solver.

### Local benchmarks

Synthetic 3-trait fixture, macOS laptop, 2026-05-04:

| Benchmark | lavaan | rust-backed | Speedup |
| --- | ---: | ---: | ---: |
| direct fixed-measurement `sem()` refit | `0.015040 s` | `0.000244 s` | `61.64x` |
| end-to-end sequential `userGWAS()` smoke, 1 SNP | `0.105400 s` | `0.013300 s` | `7.92x` |

### Next packet

1. Decide whether to broaden `userGWAS_rust()` next into `Q_SNP = TRUE`, the
   unconstrained `fix_measurement = FALSE` path, or parallel worker support.
2. If the goal is still the pure isolated-backend experiment, the most
   informative next step is probably a second real user-model family rather
   than changing outer GenomicSEM orchestration.

## 2026-05-04 13:41 PDT

### Implemented

- Expanded `userGWAS_rust()` from the first default packet to the full current
  one-factor sequential DWLS slice:
  - both `fix_measurement = TRUE/FALSE`
  - both `Q_SNP = TRUE/FALSE`
- Added native `lavaan_rust()` slot reuse for the unrestricted
  `user_gwas_dwls` model family used when `fix_measurement = FALSE`.
- Kept the supported matrix explicit in the `userGWAS_rust()` source comments
  so the strict wrapper boundary is visible where the behavior is enforced.

### Current support boundary

`userGWAS_rust()` now supports simple one-factor SNP models such as:

```r
F1 =~ A + B + C
F1 ~ SNP
```

with:

- `parallel = FALSE`
- `estimation = "DWLS"`
- `TWAS = FALSE`
- either `fix_measurement = TRUE/FALSE`
- either `Q_SNP = TRUE/FALSE`

Still unsupported:

- `parallel = TRUE`
- `TWAS = TRUE`
- ML estimation
- multi-factor or more general user-model syntax

### Validation

- `lavaanrust/tests/testthat`: 38 passing tests.
- Added native reuse coverage for the unrestricted `user_gwas_dwls` base model.
- Synthetic sequential GenomicSEM smokes, 1 SNP:
  - `fix_measurement = FALSE`, `Q_SNP = FALSE`:
    - max estimate difference: `1.04e-07`
    - max sandwich SE difference: `4.89e-08`
    - max chi-square difference: `1.16e-13`
  - `fix_measurement = FALSE`, `Q_SNP = TRUE`:
    - max `Q_SNP` difference: `3.44e-08`
    - max `Q_SNP` p-value difference: `1.65e-08`
  - `fix_measurement = TRUE`, `Q_SNP = TRUE`:
    - max `Q_SNP` difference: `2.06e-09`
    - max `Q_SNP` p-value difference: `9.86e-10`

### Local benchmarks

Synthetic 3-trait fixture, macOS laptop, 2026-05-04:

| Benchmark | lavaan | rust-backed | Speedup |
| --- | ---: | ---: | ---: |
| direct unrestricted `userGWAS()` `sem()` fit | `0.021660 s` | `0.000732 s` | `29.59x` |
| end-to-end sequential, `fix_measurement = FALSE`, `Q_SNP = FALSE` | `0.062500 s` | `0.011800 s` | `5.30x` |
| end-to-end sequential, `fix_measurement = FALSE`, `Q_SNP = TRUE` | `0.057200 s` | `0.005300 s` | `10.79x` |
| end-to-end sequential, `fix_measurement = TRUE`, `Q_SNP = FALSE` | `0.100500 s` | `0.006600 s` | `15.23x` |
| end-to-end sequential, `fix_measurement = TRUE`, `Q_SNP = TRUE` | `0.103800 s` | `0.007100 s` | `14.62x` |

### Next packet

1. Add parallel worker support for the now-complete current one-factor
   `userGWAS_rust()` slice.
2. Then decide whether to broaden `usermodel_rust()` syntax or move to a
   second real user-GWAS model family.

## 2026-05-04 15:18 PDT

### Implemented

- Added explicit parallel worker plumbing for backend-swapped GWAS wrappers:
  - `commonfactorGWAS()` and `userGWAS()` now bind the worker helper before
    entering `foreach`
  - worker package exports are selected from wrapper-local settings, so
    `_rust()` variants load `lavaanrust` workers without changing the original
    serial logic
- Enabled strict native parallel execution for:
  - `commonfactorGWAS_rust()`
  - `userGWAS_rust()`
- Added `tools/bench_lavaan_rust_parallel.R`, a reproducible benchmark harness
  that:
  - generates deterministic synthetic LDSC/GWAS inputs
  - checks sequential-vs-parallel equality for the rust-backed wrappers
  - benchmarks lavaan and rust-backed wrappers over multiple core counts

### Current support boundary

`commonfactorGWAS_rust()` now supports the current one-factor DWLS model family
with either `parallel = TRUE/FALSE`.

`userGWAS_rust()` now supports simple one-factor SNP models such as:

```r
F1 =~ A + B + C
F1 ~ SNP
```

with:

- `parallel = TRUE/FALSE`
- `estimation = "DWLS"`
- `TWAS = FALSE`
- `std.lv = FALSE`
- either `fix_measurement = TRUE/FALSE`
- either `Q_SNP = TRUE/FALSE`

Still unsupported:

- `TWAS = TRUE`
- `std.lv = TRUE`
- ML estimation
- multi-factor or more general user-model syntax

### Validation

- Local 4-SNP and benchmark-harness smoke checks matched exactly between
  rust-backed sequential and parallel output on the checked numeric columns.
- Remote panda/flex 16-CPU benchmark pod, synthetic `n_snp = 1000`, `k = 5`
  fixture:
  - `commonfactorGWAS_rust()` sequential vs 2-worker parallel max absolute
    difference: `0`
  - `userGWAS_rust()` sequential vs 2-worker parallel max absolute difference:
    `0`

### Remote benchmarks

Synthetic 1,000-SNP, 5-trait fixture on a panda/flex 16-CPU pod,
`OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`,
2026-05-04:

| Workflow | Backend | sequential | 2 cores | 4 cores | 8 cores | 16 cores |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `commonfactorGWAS()` | lavaan | `167.523 s` | `85.818 s` | `45.065 s` | `24.809 s` | `18.668 s` |
| `commonfactorGWAS_rust()` | rust-backed | `6.227 s` | `4.432 s` | `2.091 s` | `1.509 s` | `1.365 s` |
| `userGWAS()` | lavaan | `62.443 s` | `62.423 s` | `28.090 s` | `17.820 s` | `14.530 s` |
| `userGWAS_rust()` | rust-backed | `5.439 s` | `4.581 s` | `2.704 s` | `2.138 s` | `2.361 s` |

Selected comparisons:

- `commonfactorGWAS_rust()` is `26.90x` faster than lavaan sequentially and
  `13.68x` faster at 16 workers.
- `userGWAS_rust()` is `11.48x` faster than lavaan sequentially and `6.15x`
  faster at 16 workers.
- The common-factor rust path continues to improve through 16 workers on this
  fixture; the user-GWAS rust path is fastest at 8 workers here, so worker
  overhead is already visible once the per-SNP fit cost is small enough.

### Next packet

1. Broaden `usermodel_rust()` syntax or move to a second real user-GWAS model
   family.
2. If parallel performance remains a priority later, revisit task granularity
   and batching for the very fast rust-backed loops rather than only adding
   more workers.

## 2026-05-04 16:05 PDT

### Implemented

- Expanded `tools/bench_lavaan_rust_parallel.R` from the default
  fixed-measurement / `Q_SNP = TRUE` benchmark into the full current supported
  `userGWAS_rust()` matrix:
  - `fix_measurement = TRUE/FALSE`
  - `Q_SNP = TRUE/FALSE`
- Added full-matrix equivalence checks before timing:
  - lavaan vs rust-backed sequential output for every supported user-GWAS cell
  - rust-backed sequential vs 2-worker parallel output for every supported
    user-GWAS cell

### Validation

Remote panda/flex 16-CPU pod, synthetic `n_snp = 1000`, `k = 5` fixture:

| Workflow | `fix_measurement` | `Q_SNP` | Comparison | Max absolute difference |
| --- | --- | --- | --- | ---: |
| `commonfactorGWAS()` | `NA` | `NA` | rust sequential vs parallel | `0` |
| `userGWAS()` | `TRUE` | `TRUE` | lavaan vs rust sequential | `1.18e-07` |
| `userGWAS()` | `TRUE` | `TRUE` | rust sequential vs parallel | `0` |
| `userGWAS()` | `FALSE` | `TRUE` | lavaan vs rust sequential | `1.38e-07` |
| `userGWAS()` | `FALSE` | `TRUE` | rust sequential vs parallel | `0` |
| `userGWAS()` | `TRUE` | `FALSE` | lavaan vs rust sequential | `1.18e-07` |
| `userGWAS()` | `TRUE` | `FALSE` | rust sequential vs parallel | `0` |
| `userGWAS()` | `FALSE` | `FALSE` | lavaan vs rust sequential | `1.38e-07` |
| `userGWAS()` | `FALSE` | `FALSE` | rust sequential vs parallel | `0` |

### Remote benchmarks

Synthetic 1,000-SNP, 5-trait fixture on the same panda/flex 16-CPU pod,
`OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`,
2026-05-04:

| Workflow | Backend | `fix_measurement` | `Q_SNP` | sequential | 2 cores | 4 cores | 8 cores | 16 cores |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| `commonfactorGWAS()` | lavaan | `NA` | `NA` | `167.093 s` | `85.649 s` | `44.774 s` | `24.852 s` | `18.423 s` |
| `commonfactorGWAS_rust()` | rust-backed | `NA` | `NA` | `6.250 s` | `4.941 s` | `2.227 s` | `1.504 s` | `1.463 s` |
| `userGWAS()` | lavaan | `TRUE` | `TRUE` | `62.121 s` | `62.418 s` | `28.771 s` | `17.526 s` | `14.421 s` |
| `userGWAS_rust()` | rust-backed | `TRUE` | `TRUE` | `5.411 s` | `5.484 s` | `2.819 s` | `2.176 s` | `2.283 s` |
| `userGWAS()` | lavaan | `FALSE` | `TRUE` | `87.534 s` | `77.286 s` | `35.251 s` | `20.281 s` | `14.901 s` |
| `userGWAS_rust()` | rust-backed | `FALSE` | `TRUE` | `6.589 s` | `6.389 s` | `2.928 s` | `2.438 s` | `2.283 s` |
| `userGWAS()` | lavaan | `TRUE` | `FALSE` | `61.323 s` | `58.999 s` | `28.431 s` | `18.290 s` | `13.912 s` |
| `userGWAS_rust()` | rust-backed | `TRUE` | `FALSE` | `4.472 s` | `4.008 s` | `2.455 s` | `2.122 s` | `1.937 s` |
| `userGWAS()` | lavaan | `FALSE` | `FALSE` | `87.180 s` | `81.775 s` | `35.065 s` | `20.161 s` | `16.304 s` |
| `userGWAS_rust()` | rust-backed | `FALSE` | `FALSE` | `5.341 s` | `4.045 s` | `2.651 s` | `2.269 s` | `2.036 s` |

Selected comparisons:

- Sequential speedups span `11.48x` to `16.32x` across the full supported
  `userGWAS_rust()` matrix.
- At 16 workers, speedups span `6.32x` to `8.01x` across that matrix.
- `Q_SNP = FALSE` lowers rust-backed runtime materially in both measurement
  modes, especially for `fix_measurement = TRUE`.
- The best rust-backed worker count is not constant:
  - `fix_measurement = TRUE`, `Q_SNP = TRUE`: best at 8 workers
  - `fix_measurement = FALSE`, `Q_SNP = TRUE`: best at 16 workers
  - `fix_measurement = TRUE`, `Q_SNP = FALSE`: best at 16 workers
  - `fix_measurement = FALSE`, `Q_SNP = FALSE`: best at 16 workers

### Next packet

1. Broaden `usermodel_rust()` syntax or move to a second real user-GWAS model
   family.
2. If more parallel work is desired later, target batching/task granularity for
   the already-fast rust loops rather than expecting uniform gains from simply
   raising `cores`.

## 2026-05-04 16:13 PDT

### Tightening

- Added an explicit `std.lv = TRUE` rejection in `userGWAS_rust()`.
- This was already outside the intended native slice, but it was not previously
  guarded:
  - `fix_measurement = TRUE, std.lv = TRUE` failed mid-run
  - `fix_measurement = FALSE, std.lv = TRUE` completed but was not equivalent to
    lavaan
- The wrapper now keeps the experiment's strict contract: unsupported paths
  fail immediately instead of producing mixed or incorrect results.

## 2026-05-04 17:41 PDT

### Generic compiler scaffold

- Added `LAVAAN_FAST_PLAN.md` to make the next phase explicit: `lavaan_fast`
  should begin as a generic compiler layer rather than as another family of
  hand-written fast paths.
- Added a parameter-table compiler in `lavaanrust/R/compiler.R` that lowers the
  current supported subset (`=~`, `~`, `~~`) into a RAM representation:
  - directed matrix `A`
  - residual covariance matrix `S`
  - observed-variable selector `F`
- Added generic reconstruction of:
  - implied covariance `Sigma = F (I - A)^-1 S (I - A)^-T F^T`
  - analytic `d vech(Sigma) / d theta` Jacobians
- Kept unsupported operators strict-errors at the compiler boundary; `:=` still
  does not silently pass through.

### Validation

- Added compiler tests for both current user-GWAS shapes:
  - unrestricted one-factor model
  - fixed-measurement one-factor model
- The compiler now exactly reproduces the specialized rust-backed fit objects
  for both:
  - implied covariance matrices
  - Jacobian matrices exposed through `lavInspect_rust(..., "delta")`
- Local package validation after the compiler packet:
  - `R CMD INSTALL lavaanrust`
  - `Rscript -e 'testthat::test_dir("lavaanrust/tests/testthat")'`
  - result: `46` passing tests

### Next packet

1. Route one existing family through the generic compiler without changing its
   public wrapper contract.
2. Then widen the syntax/compiler surface toward the real blocker set for
   broader `userGWAS_rust()` coverage:
   - labels and equality reuse
   - direct effects
   - residual covariances
   - eventually defined parameters and nonlinear constraints

## 2026-05-04 18:50 PDT

### Native compiler evaluator

- Added `evaluate_ram_surfaces()` in Rust and the R bridge
  `.lavaan_fast_implied_surfaces_rust()`.
- The native evaluator consumes the compiler's row-wise RAM encoding and
  returns both:
  - implied observed covariance
  - analytic Jacobian over `vech(Sigma)`
- Added fixture checks proving the native evaluator matches the specialized
  rust-backed `userGWAS` surfaces for both current model shapes.

### Validation and profiling

- Local package validation:
  - `R CMD INSTALL lavaanrust`
  - `Rscript -e 'testthat::test_dir("lavaanrust/tests/testthat")'`
  - result: `50` passing tests
- Direct evaluator timing on the unrestricted user-GWAS fixture:
  - R reference evaluator, 10,000 calls: `2.274 s`
  - Rust evaluator, 10,000 calls: `1.004 s`
- I also tried routing the unrestricted `userGWAS` constructor through the
  R-side compiler surfaces directly. That was behaviorally exact but slower:
  - baseline specialized constructor path, 1,000 tiny fits: `0.842 s`
  - R compiler-backed constructor path, 1,000 tiny fits: `1.134 s`
- I did not keep that hot-path change. The next production-quality step is a
  native generic optimizer that reuses compiler-backed surfaces during fitting,
  instead of paying an extra post-fit reconstruction cost after a specialized
  optimizer has already computed them.

### Next packet

1. Add a generic native DWLS optimizer over the RAM IR.
2. Route one existing supported family through that optimizer only if the
   benchmark remains competitive with the specialized kernel.

## 2026-05-04 18:53 PDT

### Generic RAM optimizer

- Added `fit_ram_dwls()` in Rust and the R bridge
  `.lavaan_fast_fit_dwls_rust()`.
- The generic optimizer now:
  - consumes the same compiled RAM rows as the native evaluator
  - reuses generic implied covariance and Jacobian surfaces during fitting
  - clamps free diagonal covariance parameters positive
  - returns estimates, implied covariance, Jacobian, naive SEs, fit objective,
    SRMR, convergence, and iteration count
- Added a fixed-measurement `userGWAS` fixture test showing the generic optimizer
  reproduces the existing specialized rust-backed fit.

### Validation and benchmark

- Local package validation:
  - `R CMD INSTALL lavaanrust`
  - `Rscript -e 'testthat::test_dir("lavaanrust/tests/testthat")'`
  - result: `53` passing tests
- Direct low-level benchmark on the fixed-measurement user-GWAS fixture:
  - specialized kernel, 1,000 fits: `0.014 s`
  - generic RAM optimizer, 1,000 fits: `0.109 s`
- Conclusion:
  - the generic optimizer is now a viable compiler-backed coverage path
  - it is not yet competitive with the tiny analytic kernels, so those should
    remain the production hot paths for already-supported models

### Next packet

1. Broaden the compiler syntax surface enough to unlock a real currently
   unsupported model shape.
2. Use the generic optimizer there first, where the comparison is against
   unsupported behavior rather than against an already-optimal specialized
   kernel.

## 2026-05-04 18:58 PDT

### Generic parameter-table fallback

- Added a generic compiler-backed `sem_rust()` path for parameter-table models
  with free directed structure that do not match one of the existing specialized
  kernels.
- Added matching `lavaan_rust()` reuse support for the new
  `ram_dwls_generic` model kind.
- Tightened the fixed-measurement user-GWAS recognizer so extra rows no longer
  get silently accepted by the narrow specialized slice.

### Validation

- Added a parameter-table fixture with an extra free direct SNP effect
  (`A ~ SNP`) that routes through the new generic path and supports base-model
  reuse.
- Local package validation:
  - `R CMD INSTALL lavaanrust`
  - `Rscript -e 'testthat::test_dir("lavaanrust/tests/testthat")'`
  - result: `57` passing tests

### Significance

- `lavaan_fast` is no longer only a compiler experiment; it can now execute
  genuinely broader parameter-table models natively without adding a new
  handwritten solver family.
- The next limiting step for end-to-end broader `userGWAS_rust()` support is not
  the compiled parameter-table layer anymore. It is the front-end parser/model
  construction step for richer model strings before those tables exist.

## 2026-05-04 19:12 PDT

### Rank-1 Jacobian optimization

- Replaced the dense generic Jacobian path with algebraically equivalent
  outer-product updates:
  - directed-path derivatives now use columns of `B = (I - A)^-1` and rows of
    the full implied covariance
  - covariance derivatives now use outer products of columns of `B`
- Preserved equality-reuse semantics by summing contributions from every row
  sharing a free-parameter id.
- Added an explicit regression test where one free parameter controls two
  distinct directed edges.

### Validation

- Local package validation:
  - `R CMD INSTALL lavaanrust`
  - `Rscript -e 'testthat::test_dir("lavaanrust/tests/testthat")'`
  - result: `58` passing tests

### Benchmarks

Raw native evaluator, same synthetic RAM shapes before vs after:

| Observed vars | Free params | Dense Jacobian | Rank-1 Jacobian | Speedup |
| ---: | ---: | ---: | ---: | ---: |
| 4 | 10 | `9.70 us` | `3.30 us` | `2.94x` |
| 8 | 18 | `36.90 us` | `6.20 us` | `5.95x` |
| 16 | 34 | `167.67 us` | `24.00 us` | `6.99x` |
| 32 | 66 | `1330.00 us` | `132.00 us` | `10.08x` |

At fixed `8` observed variables, scaling across free-parameter counts also
improved materially:

| Free params | Dense Jacobian | Rank-1 Jacobian | Speedup |
| ---: | ---: | ---: | ---: |
| 17 | `35.20 us` | `6.07 us` | `5.80x` |
| 25 | `50.17 us` | `7.08 us` | `7.09x` |
| 41 | `79.00 us` | `9.12 us` | `8.66x` |
| 73 | `139.50 us` | `13.25 us` | `10.53x` |

The public R wrapper path improved more modestly on the tiny user-GWAS fixture:

- before: `1.004 s` for 10,000 calls
- after: `0.810 s` for 10,000 calls

`Rprof` now shows that small-model wrapper overhead is the next clear target:
matrix reconstruction, stat-name generation, and free-label generation consume
more time than the native evaluator itself.

### Next packet

1. Move more compile products into the compiled object:
   - free labels
   - stat names
   - native row groups / edge metadata
2. Add an internal flat-array path for repeated native use and wrap named R
   matrices only at the public boundary.
3. After that, revisit `fit_ram_dwls()` with an implied-only line-search path;
   the generic optimizer still pays more overhead than the specialized kernels
   on tiny models.

## 2026-05-04 20:15 PDT

### Cached compile metadata and flat surfaces

- Moved repeated compiler products into the compiled model object:
  - native free-index and fixed-value vectors
  - observed-variable indices
  - free-parameter row groups encoded as offsets plus flattened row indices
  - free-parameter labels, statistic names, and reusable dimensions
- Changed the native RAM evaluator and generic DWLS optimizer to consume the
  precompiled row-group metadata instead of rebuilding free-row groupings on
  each call.
- Added `.lavaan_fast_implied_surfaces_rust_flat()` as the internal raw-array
  path; the public helper now only adds named R matrices at the boundary.
- Added regression coverage for the cached compiler fields and for agreement
  between the flat and wrapped surface APIs.

### Validation

- Local package validation:
  - `R CMD INSTALL lavaanrust`
  - `Rscript -e 'testthat::test_dir("lavaanrust/tests/testthat")'`
  - result: `64` passing tests

### Benchmarks

Tiny unrestricted `userGWAS` fixture, 10,000 evaluator calls, local Apple
Silicon run, median of seven repetitions:

| Path | Before this packet | After this packet | Speedup |
| --- | ---: | ---: | ---: |
| Named public wrapper | `0.810 s` | `0.141 s` | `5.74x` |
| Internal flat native path | n/a | `0.079 s` | n/a |

Relative to the original generic R evaluator checkpoint from `18:50 PDT`, the
current named wrapper is now `16.13x` faster (`2.274 s -> 0.141 s`) while still
returning the same named matrices.

`Rprof` after the change shows the next small-model overhead clearly:

- wrapped path:
  - native evaluator: `51.30%` self time
  - R `matrix()` reconstruction: `28.10%`
  - list element access: `13.37%`
- flat path:
  - native evaluator: `76.54%` self time
  - list element access: `16.36%`

### Next packet

1. Revisit `fit_ram_dwls()` with an implied-only line-search path so candidate
   steps stop paying for Jacobians they do not use.
2. If repeated surface calls remain important after that, consider a reusable
   native compiled plan or external pointer before trying to shave the public
   matrix-wrapping boundary further.

## 2026-05-05 10:33 PDT

### Implied-only line-search candidates

- Split candidate-step evaluation in `fit_ram_dwls()` away from the full
  surface path.
- Candidate steps now compute only the observed implied covariance:
  - no Jacobian construction
  - no full implied covariance materialization before slicing down to observed
    variables
- Kept full implied covariance plus Jacobian construction on the current iterate
  and final fit result, where those values are actually consumed.

### Validation

- Local package validation:
  - `R CMD INSTALL lavaanrust`
  - `Rscript -e 'testthat::test_dir("lavaanrust/tests/testthat")'`
  - result: `64` passing tests

### Benchmarks

Fixed-measurement `userGWAS` generic fit fixture, 1,000 fits, local Apple
Silicon run, median of seven repetitions:

| Path | Before | After |
| --- | ---: | ---: |
| Generic RAM DWLS fit | `0.024 s` | `0.024 s` |

Synthetic unrestricted covariance fits, same harness before vs after:

| Observed vars | Fits | Before | After | Speedup |
| ---: | ---: | ---: | ---: | ---: |
| 8 | 500 | `0.053 s` | `0.049 s` | `1.08x` |
| 16 | 100 | `0.239 s` | `0.232 s` | `1.03x` |
| 24 | 30 | `0.733 s` | `0.711 s` | `1.03x` |

This is a correctness-preserving cleanup and a small win, not a major remaining
hotspot. After the previous metadata-caching packet, the generic fitter is
already close enough on tiny models that line-search candidates are no longer
the dominant cost.

### Next packet

1. Profile the generic fitter again on broader model families before adding more
   optimizer-specific micro-optimizations.
2. The next likely higher-leverage work is coverage, not another local fitter
   tweak:
   - parser/model-string support for a broader syntax subset
   - or a reusable native compiled plan if repeated surface evaluation proves
     material in a real workflow rather than only in microbenchmarks

## 2026-05-06 14:05 PDT

### Generic model-string compiler

- Added a generic parser that lowers explicit lavaan-style RAM strings directly
  into the existing parameter-table compiler path.
- Supported string syntax now includes:
  - `=~`, `~`, and `~~`
  - multiple right-hand-side terms separated by `+`
  - fixed numeric coefficients
  - `NA*` free markers
  - `start(value)*` modifiers
  - labels with equality reuse across rows
  - repeated same-edge modifier declarations such as
    `NA*A + start(.1)*A`
  - simple labeled lower bounds such as `rvA > .001`
- Kept the first generic string path deliberately explicit:
  - every observed and latent variable must still have an explicit diagonal
    variance row
  - shorthand models that rely on lavaan `auto.*` expansion remain unsupported
    until that behavior is implemented intentionally
- Generic string models now route through `sem_rust()` into the native RAM
  fitter without first requiring a hand-built parameter table.

### Validation

- Added parser regressions for modifiers, repeated terms, lower bounds, and
  generic direct-effect string fitting.
- Added a one-parameter constrained fit showing that `A ~~ rvA*A` plus
  `rvA > 1.5` is enforced by the native optimizer.
- Local package validation:
  - `R CMD INSTALL lavaanrust`
  - `Rscript -e 'testthat::test_dir("lavaanrust/tests/testthat")'`
  - result: `74` passing tests

### Coverage significance

This is the first point where broader user-authored model strings can reach the
generic Rust backend directly rather than only through handwritten specialized
parsers or preconstructed parameter tables.

The largest remaining gaps are now clearer:

1. lavaan-style automatic model expansion for omitted residual/latent variances
   and generic `std.lv` semantics
2. user-defined parameters via `:=`
3. richer symbolic constraints beyond repeated labels and simple lower bounds
4. everything outside the current DWLS covariance/RAM subset

For GenomicSEM-oriented coverage, items 1-3 are the important remaining front-end
work. Full lavaan parity is substantially larger than that and would include
means/intercepts, additional estimators, categorical machinery, groups, and much
more of lavaan's option surface.

## 2026-05-06 14:30 PDT

### Generic auto expansion and `std.lv`

- Extended the generic RAM string parser to reproduce the main lavaan auto-row
  behavior needed by shorthand measurement models:
  - observed residual variances for endogenous indicators
  - observed exogenous variances and pairwise covariances
  - latent variances
  - exogenous latent covariances
- Added generic identification behavior:
  - marker scaling when `std.lv = FALSE`
  - fixed latent variances with free loadings when `std.lv = TRUE`
- Used live lavaan `parTable()` output as the reference while implementing the
  row construction and start-value conventions, including the `0.05` latent
  variance starts on marker-scaled factors.
- Replaced the generic RAM fit object's placeholder latent-correlation matrix
  with a correlation matrix computed from the fitted latent covariance implied
  by the compiled RAM system.

### Validation

- Added regressions for:
  - shorthand auto-expansion under marker scaling
  - shorthand auto-identification under `std.lv = TRUE`
  - auto-generated exogenous observed variances/covariances
  - end-to-end fitting of a two-factor shorthand string through the generic
    native RAM path
- Local package validation:
  - `R CMD INSTALL lavaanrust`
  - `Rscript -e 'testthat::test_dir("lavaanrust/tests/testthat")'`
  - result: `88` passing tests

### Coverage significance

The generic path now covers the common shorthand measurement-model surface that
GenomicSEM users actually write, instead of only accepting fully expanded RAM
strings. The largest remaining front-end gaps are now:

1. user-defined parameters via `:=`
2. richer symbolic constraints beyond repeated labels and simple lower bounds
3. broader standardized-output compatibility for generic multi-factor models

Beyond those, reaching full lavaan coverage would still require a much larger
surface: means/intercepts, additional estimators, categorical machinery,
multi-group behavior, and other option families outside the current DWLS
covariance subset.

## 2026-05-06 14:40 PDT

### Generic defined parameters

- Added generic `:=` support on top of the RAM compiler without changing the
  native solve itself:
  - structural rows still compile into the RAM system
  - defined rows stay symbolic and are evaluated from labeled fitted parameters
  - supported expressions currently cover `+`, `-`, `*`, `/`, `^`, `sqrt()`,
    `exp()`, and `log()`
- Added a minimal `lavaan_rust_model` compatibility object for generic fits so
  unchanged GenomicSEM code can keep using `fit@Model@def.function`.
- Implemented:
  - `lav_model_get_parameters_rust()`
  - `lav_func_jacobian_complex_rust()`
- Extended generic `standardizedSolution_rust()` so defined parameters are
  re-evaluated from standardized labeled estimates.
- Tolerated row-only named covariance matrices during generic parsing. This
  matches a real `usermodel()` standardized-fit intermediate where lavaan works
  with row names even though column names have been dropped.

### Validation

- Added regressions for:
  - nested defined-parameter evaluation and model reuse
  - complex-step Jacobians used by GenomicSEM's delta-method SE code
  - unsupported defined-expression rejection
  - row-only covariance names
  - generic standardized output with re-evaluated defined parameters
- Local package validation:
  - `R CMD INSTALL lavaanrust`
  - `Rscript -e 'testthat::test_dir("lavaanrust/tests/testthat")'`
  - result: `98` passing tests
- End-to-end smoke:
  - `usermodel_rust()` now completes for a labeled one-factor model with
    `omega := l2 * l3`
  - the output contains unstandardized, `STD_Genotype`, and `STD_All` values for
    the defined parameter.
  - max absolute old-vs-rust difference across those defined-parameter outputs:
    `5.78e-08`

### Coverage significance

This closes the largest remaining front-end gap for ordinary user-authored DWLS
RAM models. The next meaningful symbolic gap is richer constraint syntax, not
basic model-string coverage.

## 2026-05-06 16:05 PDT

### userGWAS std.lv coverage

- Re-opened `userGWAS_rust()` for `std.lv = TRUE` instead of keeping the earlier
  defensive rejection.
- Kept the specialized one-factor SNP solver on the `std.lv = FALSE` slice and
  route `std.lv = TRUE` through the generic RAM evaluator, where the broader
  syntax support already exists.
- Matched lavaan's parameter-table behavior in three places that matter for the
  fixed-measurement path:
  - auto-generated row order for endogenous residuals, latent variances, and
    exogenous observed rows
  - a `lower` column on specialized `parTable_rust()` output
  - lavaan-style renumbering of unlabeled free rows when a parameter table is
    re-ingested after GenomicSEM reconstructs the with-SNP model
- Added top-level `GenomicSEM` regression coverage for:
  - the full `userGWAS_rust()` matrix over `std.lv = TRUE/FALSE`,
    `fix_measurement = TRUE/FALSE`, and `Q_SNP = TRUE/FALSE`
  - a broader end-to-end model using labels, a direct SNP effect, and
    `combo := gamma * direct`

### Validation

- Local regression suites:
  - `lavaanrust/tests/testthat`: `104` passing tests
  - top-level `GenomicSEM` wrapper tests: `22` passing expectations
- Simple one-factor `userGWAS()` equivalence fixture, 2 SNPs:
  - all 8 supported matrix cells preserve row order exactly
  - maximum old-vs-rust absolute difference across the matrix: `1.00e-07`
- Flexible direct-effect fixture, 2 SNPs, `std.lv = TRUE`:
  - `fix_measurement = FALSE`: maximum absolute difference `1.10e-08`
  - `fix_measurement = TRUE`: maximum absolute difference `1.19e-08`

### Coverage significance

`userGWAS_rust()` now covers the full current simple one-factor DWLS matrix with
either loading-identification convention, plus a meaningfully broader generic
workflow with labeled parameters and a defined parameter. That makes the
generic compiler path part of the real wrapper contract rather than just a
standalone evaluator experiment.

## 2026-05-06 17:06 PDT

### Broader surface coverage

- Promoted the ad hoc multi-factor probes into top-level wrapper regressions:
  - `userGWAS_rust()` now has formal two-factor coverage across
    `std.lv = TRUE/FALSE`, `fix_measurement = TRUE/FALSE`, and
    `Q_SNP = TRUE/FALSE`
  - `usermodel_rust()` now has formal two-factor coverage under both
    identification conventions
- Added `tools/bench_lavaan_rust_surface_matrix.R`, a sequential benchmark
  harness that covers:
  - one-factor and two-factor `usermodel()` cases
  - the full simple one-factor `userGWAS()` matrix
  - flexible one-factor `userGWAS()` with labels, a direct SNP effect, and `:=`
  - the full two-factor `userGWAS()` matrix

### Richer constraints

- Expanded the generic parser from lower-bound-only support to simple symbolic
  box constraints and label equalities:
  - lower bounds via `>`
  - upper bounds via `<`
  - explicit equalities such as `l2 == l3`
- Generalized the Rust RAM optimizer from lower-bound projection to box-bound
  projection.
- Kept the boundary deliberate: this packet supports scalar bounds and simple
  label equalities, not arbitrary nonlinear constraint expressions.

### Validation

- Local regression suites:
  - `lavaanrust/tests/testthat`: `116` passing tests
  - top-level `GenomicSEM` wrapper tests: `42` passing expectations
- Two-factor end-to-end `userGWAS()` smoke before formalizing the tests:
  - all 8 option cells preserved row order exactly
  - max old-vs-rust absolute difference across the matrix: `6.12e-08`
- New constraint regressions cover:
  - box-bound parsing and propagation across equality-linked labels
  - upper-bound enforcement in the generic fitter
  - rejection of incompatible bound sets

### Local benchmark smoke

`tools/bench_lavaan_rust_surface_matrix.R n_snp=20 repeats=1`,
macOS laptop, 2026-05-06:

| Workflow surface | lavaan | rust-backed | Speedup |
| --- | ---: | ---: | ---: |
| one-factor `usermodel()`, `std.lv = FALSE` | `0.064 s` | `0.005 s` | `12.80x` |
| two-factor `usermodel()`, `std.lv = FALSE` | `0.072 s` | `0.022 s` | `3.27x` |
| simple one-factor `userGWAS()`, `std.lv = FALSE`, `fix_measurement = FALSE`, `Q_SNP = FALSE` | `0.433 s` | `0.040 s` | `10.83x` |
| flexible one-factor `userGWAS()`, `std.lv = TRUE`, `fix_measurement = FALSE` | `0.406 s` | `0.058 s` | `7.00x` |
| two-factor `userGWAS()`, `std.lv = FALSE`, `fix_measurement = FALSE`, `Q_SNP = FALSE` | `0.508 s` | `0.068 s` | `7.47x` |

This is a surface-coverage smoke, not the release benchmark suite. Its purpose
is to keep broader supported paths measured while the remote scaling benchmark
continues to cover the heavier one-factor GWAS workloads.

### Coverage significance

The branch is no longer only a one-factor backend experiment. The generic path
now has explicit wrapper-level coverage for a broader multi-factor workflow and
supports the next useful slice of lavaan's symbolic layer while still keeping
unsupported nonlinear constraints out of scope.

## 2026-05-07 10:50 PDT

### Constrained wrapper coverage

- Added bounded-model regressions through the public wrapper layer:
  - `usermodel_rust()` now has end-to-end coverage for a labeled residual
    variance lower bound.
  - `userGWAS_rust()` now has end-to-end coverage for a labeled loading lower
    bound.
- Normalized the wrapper tests so they compare the modeled structural result
  rather than incidental row-name representation or whether a given installed
  lavaan build echoes textual bound rows back into the output table.
- Kept repeated-label equality constraints out of the wrapper regression for
  this packet after profiling exposed a real compatibility gap:
  - the generic compiler can fit them
  - but the public wrappers still do not fully reproduce lavaan's exposed
    free-parameter / equality-row shape for those models
  - that should be fixed explicitly rather than hidden inside a loose test

### Broader workflow profiling

Added `tools/profile_lavaan_rust_workflows.R` so broader supported wrappers can
be profiled one case at a time after a warm-up run:

- `usermodel_bounded`
- `usergwas_flexible_one_factor`
- `usergwas_two_factor`

Remote panda/flex 16-CPU profiles on the 20-SNP fixture show the next clear
optimization target:

- the wrapper-level `userGWAS` machinery still dominates total time
- the Rust fitter itself is only a small slice of the broad generic wrapper
  profiles on this small fixture
- parsing and data-frame assembly remain visible in both the bounded
  `usermodel()` case and the flexible `userGWAS()` case

That means the next high-impact work should be plan reuse / wrapper assembly,
not another local optimizer micro-tweak.

### Remote surface smoke

Remote panda/flex 16-CPU pod, `tools/bench_lavaan_rust_surface_matrix.R
n_snp=20 repeats=1`:

| Workflow | Representative case | lavaan | lavaanrust | Speedup |
| --- | --- | ---: | ---: | ---: |
| `usermodel()` | one-factor, `std.lv = FALSE` | `0.140 s` | `0.014 s` | `10.00x` |
| `usermodel()` | two-factor, `std.lv = FALSE` | `0.164 s` | `0.059 s` | `2.78x` |
| `userGWAS()` | simple one-factor, `fix_measurement = FALSE`, `Q_SNP = FALSE` | `1.521 s` | `0.114 s` | `13.34x` |
| `userGWAS()` | flexible one-factor, `fix_measurement = FALSE` | `1.459 s` | `0.158 s` | `9.23x` |
| `userGWAS()` | two-factor, `fix_measurement = FALSE`, `Q_SNP = FALSE` | `1.725 s` | `0.178 s` | `9.69x` |

The surface smoke remains numerically close to lavaan across the broadened
fixtures; the largest reported absolute difference in that run was
`3.09e-05`.

### Parameter-count scaling

Added `tools/bench_lavaan_rust_parameter_scaling.R`, which measures two related
things on a widening full-covariance family:

1. end-to-end fit time for lavaan vs lavaanrust
2. raw Rust implied-surface evaluation time with the compiled model held fixed

The widening family spans `10` to `528` free parameters. On the remote
panda/flex 16-CPU pod:

| Surface | Empirical exponent from `lm(log(time) ~ log(n_free))` | `R^2` |
| --- | ---: | ---: |
| lavaan full fit | `0.51` | `0.875` |
| lavaanrust full fit | `1.08` | `0.994` |
| lavaanrust raw surface evaluator | `1.39` | `0.939` |

Selected median timings from the same run:

| Free params | lavaan full fit | lavaanrust full fit | lavaanrust surface eval |
| ---: | ---: | ---: | ---: |
| `10` | `0.026 s` | `0.014 s` | `0.000016 s` |
| `78` | `0.040 s` | `0.101 s` | `0.000074 s` |
| `210` | `0.075 s` | `0.299 s` | `0.000460 s` |
| `528` | `0.195 s` | `1.051 s` | `0.002862 s` |

This is an empirical widening-family study, not a symbolic proof of asymptotic
complexity. The useful conclusion is narrower: over the range relevant to the
current generic evaluator, the Rust backend is not hiding a sudden cubic
parameter-count cliff, while the raw surface path is already the place where
superlinear growth becomes visible as model width increases.

### Matched nonsaturated scaling follow-up

Added `tools/bench_lavaan_rust_matched_scaling.R` to answer the more specific
question the full-covariance benchmark cannot answer cleanly: how do lavaan and
lavaanrust scale on the same nonsaturated generic family when variable count and
free-parameter count are varied separately?

The benchmark has two families:

1. `variable_scaling`: a balanced two-factor CFA with no residual covariances,
   widened from `6` to `48` observed variables. This grows the observed
   covariance dimension while keeping the model nonsaturated and staying on the
   Rust generic RAM path.
2. `parameter_scaling`: the same two-factor CFA with `24` observed variables
   held fixed while residual covariance parameters are added from `0` to `176`.
   This grows free-parameter count without changing the observed dimension or
   the data matrix.

Both backends fit the exact same model/data pair at each point, and the script
records the Rust model kind so accidental dispatch into a specialized fast path
is visible rather than inferred after the fact.

Remote panda/flex 16-CPU run, 2026-05-07:

| Family | Axis | lavaan exponent | lavaanrust exponent | lavaanrust raw surface exponent |
| --- | --- | ---: | ---: | ---: |
| `variable_scaling` | `n_vars` | `0.73` | `1.37` | `1.87` |
| `variable_scaling` | `n_free` | `0.76` | `1.42` | `1.94` |
| `parameter_scaling` | `n_free` | `0.66` | `1.14` | `0.90` |

Selected medians from the same run:

| Family | Shape | lavaan full fit | lavaanrust full fit | lavaanrust raw surface |
| --- | --- | ---: | ---: | ---: |
| `variable_scaling` | `6` vars, `13` free | `0.033 s` | `0.017 s` | `0.000018 s` |
| `variable_scaling` | `24` vars, `49` free | `0.056 s` | `0.059 s` | `0.000125 s` |
| `variable_scaling` | `48` vars, `97` free | `0.154 s` | `0.265 s` | `0.000821 s` |
| `parameter_scaling` | `24` vars, `49` free | `0.055 s` | `0.059 s` | `0.000125 s` |
| `parameter_scaling` | `24` vars, `121` free | `0.082 s` | `0.184 s` | `0.000269 s` |
| `parameter_scaling` | `24` vars, `225` free | `0.154 s` | `0.364 s` | `0.000475 s` |

All rows in that run stayed on `ram_dwls_generic`, converged in both backends,
and matched fitted covariance outputs to at most `2.54e-07`.

The matched benchmark changes the interpretation of the earlier saturated
family sweep. On the generic nonsaturated path, there is not yet evidence that
lavaanrust has a better asymptotic class than lavaan; if anything, its
end-to-end full-fit curve is steeper over the current range. The large practical
speedups elsewhere in the branch therefore come from specialized kernels,
batched orchestration, lower interpreter overhead, and parallel work, not from
the generic RAM fallback already having superior scaling.

## 2026-05-07 11:18 PDT

### Matched-scaling diagnosis and fixes

The matched benchmark exposed two avoidable costs in the generic path:

1. The Rust DWLS fitters treated `WLS.V` as dense even when it was the diagonal
   matrix supplied by the DWLS workflows in these benchmarks. That made each
   iteration pay dense `W * delta` work.
2. The generic R parser built one-row `data.frame`s inside a loop and then
   row-bound them. `Rprof` on the `48`-variable case showed `data.frame()` and
   its conversion helpers dominating parser time.

Changes made:

- `b47c47c` adds a shared `DwlsWeights` representation in Rust. Diagonal
  matrices now use elementwise row scaling while genuinely dense matrices keep
  the general path.
- `990a133` builds the parser's structural parameter table column-wise and
  creates one `data.frame` at the end.

Remote panda/flex component profile on the noisy `48`-variable matched model:

| Stage | Before parser fix | After parser fix |
| --- | ---: | ---: |
| generic parser | `0.091 s` | `0.009 s` |
| Rust fit core | `0.029 s` | `0.029 s` |
| end-to-end `sem_rust()` | `0.140 s` | `0.050 s` |

Remote matched-scaling progression:

| State | `48` vars / `97` free | variable-family exponent vs `n_free` | `24` vars / `225` free | parameter-family exponent vs `n_free` |
| --- | ---: | ---: | ---: | ---: |
| before | `0.265 s` | `1.42` | `0.364 s` | `1.14` |
| after diagonal DWLS weights | `0.129 s` | `1.03` | `0.318 s` | `1.12` |
| after parser table vectorization | `0.051 s` | `1.13` | `0.128 s` | `1.51` |

The last exponent increase in the parameter family is not a regression in wall
time. Removing large fixed parser overhead exposes the real remaining generic
solver cost more clearly: with `n_vars` fixed, the raw surface evaluator grows
roughly linearly in free parameters, while the Rust fit core grows much faster
because every Gauss-Newton iteration still forms and solves dense normal
equations. In the remote parameter-family component profile, `fit_core` grows
from `0.0025 s` at `49` free parameters to `0.0603 s` at `225`.

That is now the next substantive generic-path question. The two incidental
scaling penalties are fixed; further improvement requires changing the optimizer
strategy or exploiting more structure than the current fully generic dense
normal-equation path does.

### Validation

- Local:
  - `Rscript -e 'pkgload::load_all(quiet=TRUE); testthat::test_dir("lavaanrust/tests/testthat"); testthat::test_dir("tests/testthat")'`
  - result: `116` lavaanrust tests and `46` wrapper tests passing
- Remote panda/flex 16-CPU pod:
  - `Rscript -e 'library(GenomicSEM); testthat::test_dir("lavaanrust/tests/testthat"); testthat::test_dir("tests/testthat")'`
  - result: `116` lavaanrust tests and `46` wrapper tests passing

## 2026-05-07 11:30 PDT

### Sparse RAM Jacobian accumulation

The next generic-path change targets the RAM surface evaluator itself rather than
wrapper or parser overhead:

- `c823935` rewrites generic RAM Jacobian assembly around the rank-one structure
  of each RAM row derivative.
- The evaluator now carries only the observed rows of `(I - A)^-1`, computes the
  observed-by-variable covariance projection once, and accumulates lower-triangle
  Jacobian entries through sparse outer-product support sets.
- This keeps the evaluator generic: dense models still take the dense work, but
  exact zeros in the RAM projection are no longer revisited for every free
  parameter/statistic pair.
- Added Rust unit tests that compare sparse lower-triangle accumulation against
  the dense symmetric outer-product formulas.

Remote matched-scaling progression on the same panda/flex 16-CPU pod:

| State | `48` vars / `97` free full fit | raw surface | variable-family surface exponent vs `n_free` | `24` vars / `225` free full fit | raw surface | parameter-family surface exponent vs `n_free` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| after parser table vectorization | `0.051 s` | `0.000801 s` | `1.923` | `0.128 s` | `0.000466 s` | `0.867` |
| after sparse RAM Jacobians | `0.040 s` | `0.000237 s` | `1.273` | `0.124 s` | `0.000114 s` | `0.429` |

That is a `3.4x` raw-surface win on the variable-growth family and a `4.1x`
raw-surface win on the free-parameter-growth family. The end-to-end generic fit
improves too, but less dramatically: `1.28x` on the `48`-variable case and only
`1.03x` on the `225`-free-parameter case. The residual parameter-family wall is
therefore still outside Jacobian construction: dense Gauss-Newton normal
equation formation/solve remains the next structural optimizer cost.

### Validation

- Local:
  - `cargo test --manifest-path lavaanrust/src/rust/Cargo.toml`
  - `Rscript -e 'pkgload::load_all(quiet=TRUE, recompile=TRUE); testthat::test_dir("lavaanrust/tests/testthat"); testthat::test_dir("tests/testthat")'`
  - result: `4` Rust tests, `116` lavaanrust tests, and `46` wrapper tests passing
- Remote panda/flex 16-CPU pod:
  - `R CMD INSTALL lavaanrust`
  - `Rscript -e 'library(GenomicSEM); testthat::test_dir("lavaanrust/tests/testthat"); testthat::test_dir("tests/testthat")'`
  - `Rscript tools/bench_lavaan_rust_matched_scaling.R`
  - result: `116` lavaanrust tests and `46` wrapper tests passing; benchmark rows
    still matched fitted covariance outputs to at most `2.53e-07`

## 2026-05-07 11:47 PDT

### Generic optimizer cleanup after sparse Jacobians

The sparse-Jacobian work exposed two remaining generic-fit issues:

1. `sem_rust(..., se = "none")` still paid for the final dense bread inverse,
   while the matched lavaan arm actually skipped SE work. That made the previous
   `se = "none"` comparison slightly unfair to the Rust path.
2. Even with a sparse Jacobian, the generic optimizer still formed `J'WJ`
   densely for diagonal-DWLS fits.

Changes made:

- propagate `se = "none"` into the generic RAM fitter and skip the unused final
  bread inversion when naive SEs are not requested
- solve damped Gauss-Newton systems with Cholesky first and keep the previous LU
  path as a fallback
- build diagonal-DWLS normal equations row-sparsely when exact Jacobian row
  support makes that materially cheaper than a dense crossproduct
- add regression coverage for `compute_se = FALSE`, sparse diagonal normal
  equations, and the Cholesky/LU solver helper

Remote panda/flex matched-scaling progression on the fixed-`24`-variable family:

| State | `225` free full fit | parameter-family exponent vs `n_free` |
| --- | ---: | ---: |
| after sparse RAM Jacobians | `0.124 s` | `1.509` |
| after honoring `se = "none"` | `0.115 s` | `1.473` |
| after sparse diagonal normal equations | `0.096 s` | `1.380` |
| after using Cholesky in the actual generic loop | `0.083 s` | `1.285` |

The last row is a repeat run on the same installed build; the first sparse-normal
run measured `0.097 s` and exponent `1.421`, so the gain is stable at the packet
level even though the pod shows normal benchmark jitter. On the widening
variable-count family, the largest Rust fit is now `0.042 s` in the repeat run.

This fixes the benchmark confounder and removes another avoidable dense step,
but it does not remove the remaining asymptotic wall. A focused instrumented run
on the `225`-free-parameter case measured `22` generic iterations and, inside
the native core, about `1.0 ms` total in surface construction, `2.0 ms` in
normal-equation formation, `10.1 ms` in solves, and `0.4 ms` in candidate
evaluation. I also tested a capped matrix-free conjugate-gradient solve with
dense fallback; it preserved outputs but regressed the largest case to `0.089 s`,
so it was dropped rather than carried forward. A follow-up diagonal-fringe
elimination experiment matched dense solves but left the largest case effectively
unchanged at `0.084 s`, so that extra machinery was dropped too. The generic
fallback still factors a dense damped system each iteration; that dense solve is
now the next structural optimization target if we want the free-parameter
scaling curve to flatten further.

### Validation

- Local:
  - `cargo test --manifest-path lavaanrust/src/rust/Cargo.toml`
  - `Rscript -e 'pkgload::load_all("lavaanrust", quiet=TRUE, recompile=TRUE); testthat::test_dir("lavaanrust/tests/testthat"); pkgload::load_all(quiet=TRUE); testthat::test_dir("tests/testthat")'`
  - result: `6` Rust tests, `119` lavaanrust tests, and `46` wrapper tests passing
- Remote panda/flex 16-CPU pod:
  - `R CMD INSTALL lavaanrust`
  - `Rscript -e 'library(GenomicSEM); testthat::test_dir("lavaanrust/tests/testthat"); testthat::test_dir("tests/testthat")'`
  - `Rscript tools/bench_lavaan_rust_matched_scaling.R`
  - result: `119` lavaanrust tests and `46` wrapper tests passing; matched
    benchmark rows still agree with lavaan to at most `2.53e-07`

## 2026-05-08 00:22 PDT

### Native plan reuse for repeated generic fits

The generic slot-reuse path was still reconstructing the whole compiled RAM
representation from `slotModel@par_table` on every repeated fit. That was
unnecessary work in GWAS-style loops where the model structure is fixed and only
the observed moments change.

Changes made:

- added a reusable native `RamDwlsPlan` external pointer that keeps immutable
  RAM row structure, support groups, bounds, and dimensions in Rust
- kept free starting values outside the plan so repeated fits preserve the
  existing base-model semantics
- stored both the serializable R compiled object and the native plan on generic
  `lavaan_rust_model` objects
- rebuilt the native plan lazily if an R worker serialization step invalidates
  the external pointer
- added regression coverage for ordinary plan reuse and for the serialized-plan
  rebuild path
- added `tools/bench_lavaan_rust_plan_reuse.R` to compare the old repeated
  recompile path against native plan reuse without changing fit-object assembly

Remote panda/flex 16-CPU microbenchmark, same generic model and same
fit-object reconstruction on both arms:

| Repeated fits | recompile each fit | native plan reuse | Speedup |
| ---: | ---: | ---: | ---: |
| `10` | `0.016 s` | `0.008 s` | `2.00x` |
| `100` | `0.157 s` | `0.083 s` | `1.89x` |
| `1000` | `1.571 s` | `0.818 s` | `1.92x` |

The maximum estimate difference in every measured row was exactly `0`. This is a
useful standalone repeated-fit improvement, but it does not change the larger
wrapper shape: `userGWAS_rust()` and `commonfactorGWAS_rust()` still make one
backend call and rebuild one R fit object per SNP. A future batched worker path
would be the next larger orchestration win if we decide to move beyond strict
one-SNP-at-a-time reuse.

### Validation

- Local:
  - `cargo test --manifest-path lavaanrust/src/rust/Cargo.toml`
  - `Rscript -e 'pkgload::load_all("lavaanrust", quiet=TRUE, recompile=TRUE); testthat::test_dir("lavaanrust/tests/testthat"); pkgload::load_all(quiet=TRUE); testthat::test_dir("tests/testthat")'`
  - result: `6` Rust tests, `128` lavaanrust tests, and `46` wrapper tests
    passing
- Remote panda/flex 16-CPU pod:
  - `R CMD INSTALL lavaanrust`
  - `Rscript -e 'library(GenomicSEM); testthat::test_dir("lavaanrust/tests/testthat"); testthat::test_dir("tests/testthat")'`
  - `Rscript tools/bench_lavaan_rust_plan_reuse.R counts=1,10,100,1000 repeats=5`
  - result: `128` lavaanrust tests and `46` wrapper tests passing; repeated-fit
    plan reuse stayed estimate-identical to recompiling

## 2026-05-08 00:36 PDT

### Repeated-fit scaling across lavaan and lavaanrust

Added `tools/bench_lavaan_rust_repeated_fits.R` to measure the repeated-fit
workload directly instead of inferring from single-fit scaling. The benchmark
uses the same flexible one-factor GWAS-style model family across all arms,
perturbs only the SNP-side observed moments between fits, and compares:

1. `lavaan` slot reuse, mirroring the old GenomicSEM loop
2. the pre-plan-reuse generic lavaanrust behavior, simulated by forcing the
   current build through `.fit_lavaan_fast_ram_model()` on every fit
3. the new native-plan reuse path through `lavaan_rust()`

Remote panda/flex 16-CPU result:

| Repeated fits | `lavaan` | lavaanrust recompile | lavaanrust plan reuse | Plan reuse vs lavaan | Plan reuse vs recompile |
| ---: | ---: | ---: | ---: | ---: | ---: |
| `10` | `0.329 s` | `0.018 s` | `0.011 s` | `29.91x` | `1.64x` |
| `100` | `3.623 s` | `0.179 s` | `0.105 s` | `34.50x` | `1.70x` |
| `1000` | `36.935 s` | `1.797 s` | `1.036 s` | `35.65x` | `1.73x` |
| `10000` | `370.877 s` | `18.866 s` | `10.643 s` | `34.85x` | `1.77x` |

The `1` through `1000` rows were measured with `3` repeats and summarized by
the median. The `10000` row is one direct measurement because the lavaan arm
alone takes more than six minutes there; it agrees closely with the linear
trend from the smaller repeated rows. Median slope over the `10` through `1000`
rows was:

| Backend | slope |
| --- | ---: |
| `lavaan` | `36.991 ms / fit` |
| lavaanrust recompile | `1.797 ms / fit` |
| lavaanrust plan reuse | `1.035 ms / fit` |

Correctness stayed tight throughout:

- maximum fitted-covariance difference versus lavaan was at most
  `1.073222e-09`
- the recompiling and plan-reuse lavaanrust arms had exactly identical
  estimates in every reported row

The repeated-fit benchmark answers the next optimization question cleanly.
Native plan reuse removes a real chunk of repeated setup work, but it is no
longer the dominant cost. A follow-up 10,000-fit diagnostic on the same pod
measured:

| Path | `10000` fits |
| --- | ---: |
| native plan fitter only | `0.677 s` |
| full plan-reuse `lavaan_rust()` path | `10.593 s` |

So after plan reuse, the one-fit-at-a-time R wrapper / fit-object reconstruction
path is about `15.65x` larger than the native fitter itself. That makes batching
the repeated-fit wrapper path, rather than more native fitter micro-tuning, the
next high-impact optimization if we are willing to diverge from strict
one-SNP-at-a-time GenomicSEM control flow.
