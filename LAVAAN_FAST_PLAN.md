# lavaan_fast generic compiler plan

## Why start here

The specialized `lavaanrust` solvers have already shown that replacing lavaan
under unchanged GenomicSEM wrappers can produce large gains. The next risk is
coverage: if every new model family needs its own parser, parameter-table
builder, and solver, the backend will become a collection of bespoke fast paths.

`lavaan_fast` should therefore start as a generic compiler layer. The first
compiler boundary is the lavaan-style parameter table, not raw syntax:

- GenomicSEM already consumes and mutates parameter tables.
- A parameter table carries free/fixed status and equality structure directly.
- Syntax parsing and model compilation are separable problems.
- Compiling parameter tables first lets us reuse the existing fixtures while we
  test the generic numerical representation.

## Initial representation

The first compiler packet lowers a supported parameter table into RAM form:

- `A`: directed paths, including both `~` regressions and `=~` loadings
- `S`: residual covariance matrix from `~~` rows
- `F`: selector from the full variable system to observed variables

The implied observed covariance is:

```text
Sigma = F (I - A)^-1 S (I - A)^-T F^T
```

This representation is broad enough to cover the current one-factor
`userGWAS()` families and is the natural bridge toward multi-factor path models,
direct SNP effects, residual covariances, and equality-constrained parameter
tables.

## First compiler packet

The first compiler packet is deliberately non-invasive:

1. Compile existing fitted parameter tables into a generic RAM intermediate
   representation.
2. Recompute implied covariance matrices from that IR.
3. Recompute analytic Jacobians over `vech(Sigma)` from that IR.
4. Prove exact agreement against the current rust-backed fitted objects for:
   - unrestricted one-factor `userGWAS()`
   - fixed-measurement one-factor `userGWAS()`
5. Keep unsupported operators such as `:=` strict-errors for now.

No existing fit path is switched over in this packet. The point is to establish
the compiler seam before moving solvers onto it.

## Next compiler packets

1. Use the native generic DWLS optimizer as the coverage path for model shapes
   that do not yet have a specialized kernel.
2. Keep the existing specialized kernels for hot paths where they are materially
   faster than the generic solver.
3. Keep the generic evaluator cheap to reuse repeatedly:
   - done: cache free labels, stat names, native row groups, and dimensions in
     the compiled object
   - done: add an internal flat-array surface path and wrap named R matrices
     only at the public boundary
   - done: detect diagonal DWLS weights and avoid dense weight-matrix products
     in Rust fitters
   - done: exploit sparse rank-one RAM Jacobian structure so exact-zero support
     entries are not revisited during lower-triangle accumulation
   - done: keep a reusable native RAM plan on generic fit objects, while
     retaining a serializable compiled R fallback for worker processes that
     cannot preserve external pointers
   - done: benchmark repeated generic fits through `10,000` iterations; native
     plan reuse is useful, but the remaining one-fit-at-a-time R wrapper path is
     now much larger than the native fitter itself
   - next: add repeated-fit batching if we are willing to move the `_rust()`
     wrappers away from the current one-SNP-at-a-time GenomicSEM control flow
4. Keep the generic optimizer simple and measured:
   - done: add an implied-only line-search path so candidate steps do not rebuild
     Jacobians they never consume
   - done: profile broader supported wrapper workflows before pursuing more
     local optimizer micro-optimizations
   - done: add a matched nonsaturated scaling benchmark that separates
     variable-count growth from free-parameter growth
   - done: vectorize generic parser table construction after the matched
     benchmark showed row-wise `data.frame()` assembly dominating end-to-end
     generic fits
   - done: honor `se = "none"` in the generic RAM path so benchmarks and callers
     do not pay for an unused final bread inverse
   - done: exploit sparse diagonal-DWLS normal-equation formation when exact
     Jacobian row sparsity makes it cheaper than dense crossproducts
   - done: use Cholesky for the actual generic damped solve path, with LU kept as
     a fallback
   - next: evaluate whether the generic fallback should keep forming dense
     Gauss-Newton normal equations for large free-parameter counts, or switch to
     a better solver strategy once the supported surface is broad enough; a
     first capped CG experiment was slower on the matched family and was dropped
5. Add parser support for a larger syntax subset:
   - done: explicit RAM strings using `=~`, `~`, and `~~`
   - done: fixed coefficients, `NA*`, `start()`, labels/equality reuse,
     residual covariances, direct effects, and simple labeled lower bounds
   - done: lavaan-style auto expansion for omitted residual/latent variances,
     exogenous covariances, and generic `std.lv` behavior
   - done: user-defined parameters via `:=` for the supported arithmetic subset
   - next: decide how much more of lavaan's symbolic layer belongs in the fast
     path
6. Extend the symbolic layer only where GenomicSEM needs it next:
   - done: equality reuse via repeated labels
   - done: simple lower-bound constraints
   - done: simple box bounds plus explicit label equalities
   - done: `:=` evaluation plus the lavaan-compatible free-parameter/Jacobian
     helpers GenomicSEM uses for delta-method SEs
   - next: finish lavaan-compatible wrapper exposure for repeated-label
     equality constraints before broadening into richer nonlinear constraints.

## Current parser boundary

The generic string path now accepts fully explicit RAM-style models such as:

```text
F1 =~ 1*A + l2*B + l3*C
F1 ~ beta*SNP
A ~ direct*SNP
B ~ direct*SNP
A ~~ rvA*A
B ~~ rvB*B
C ~~ rvC*C
F1 ~~ psi*F1
SNP ~~ 0.42*SNP
indirect := beta * direct
rvA > .001
```

That covers multi-factor/direct-effect model strings and the common shorthand
case where lavaan auto-generates omitted residual/latent variances and
exogenous covariances. The generic path now supports both marker-scaled and
`std.lv = TRUE` identification for shorthand measurement models, plus
user-defined parameters built from labeled parameters with `+`, `-`, `*`, `/`,
`^`, `sqrt()`, `exp()`, and `log()`.

The main remaining gaps to broad GenomicSEM workflow coverage are:

1. nonlinear or expression-valued constraint syntax beyond simple label
   equalities and scalar box bounds
2. any defined-parameter expressions outside the current arithmetic/function
   subset
3. non-DWLS estimators and the wider parts of lavaan outside the RAM covariance
   subset

The decision point after that is whether `lavaan_fast` remains a compiler that
feeds a few specialized kernels plus a generic fallback, or becomes the front
end for a fully generic SEM optimizer.
