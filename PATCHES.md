# Patches applied to upstream SQIsign-Fixed-Precision

This is a patched copy of the upstream **SQIsign-Fixed-Precision** baseline used
to produce the KLKL25 numbers in our bench. The patches below align
`RandomIdealGivenNorm` and ideal construction with the SQIsign specification,
enforce the fixed-precision assumptions, and address build and arithmetic
correctness issues.

## Precision budgets

Section 3 with the Round-3 primes gives the following magnitude budgets.
The storage is signed two's complement with 64-bit limbs.

| Level | Paper | Line-13 expression | Limbs | Raw width | Signed capacity |
| --- | ---: | ---: | ---: | ---: | ---: |
| I | 9,555 | 9,565 | 150 | 9,600 | 9,599 |
| III | 14,572 | 14,569 | 228 | 14,592 | 14,591 |
| V | 19,145 | 19,145 | 300 | 19,200 | 19,199 |

The production and test assertions cover both Section 3 expressions.
[The branch notes](docs/new-parameter.md) give the Level-I correction and the
scope of the code-level bound checks.

## Patch list

### 1. `intbig.h` — per-variant `IBZ_LIMBS` dispatch

`src/quaternion/ref/generic/include/intbig.h`

The pristine upstream uses

```c
#if SQISIGN_VARIANT == SQISIGN_LVL1
#define IBZ_LIMBS 110
#elif SQISIGN_VARIANT == SQISIGN_LVL3
#define IBZ_LIMBS 168
...
```

but `SQISIGN_VARIANT` is supplied by the build system as a bare token
(`lvl1`/`lvl3`/`lvl5`), not a defined macro. At preprocessor `#if`-time,
identifier tokens that aren't macros are substituted with 0, so every branch
collapses to `0 == 0` and the lvl1 `IBZ_LIMBS = 110` is silently picked for
every variant. The lvl3/lvl5 precomp arrays (encoded with 168/222 limbs)
then overflow the 110-limb buffer and refuse to compile, or worse, succeed
silently when array bounds aren't checked.

Fix: a two-level `_SQI_LIMBS / _SQI_LIMBS_` token-paste indirection that
forces `SQISIGN_VARIANT` to be expanded before `##`-pasting it onto
`_SQI_LIMBS_`. The result is real macro substitution, not the broken
token-equality `#if`. The current sizes above match the regenerated precomputation.

### 2. `sign.c` — silence 8 `printf("--<step> done--\n")` traces

`src/signature/ref/lvlx/sign.c`

Each call to `protocols_sign` writes 8 `--<step> done--` lines to stdout.
During a sign-retry storm these can dominate stdout and drift the measured
cycle counts by roughly 10 %. Replaced each with `(void)0;`. Behavior of
the sign protocol is unchanged.

### 3. `sign.c` — robust, silent retry cleanup

`src/signature/ref/lvlx/sign.c`

The signing retry loop no longer emits diagnostic output on `stderr`, so the
library stays silent. Auxiliary ideals are finalized on both sampler success
and failure, avoiding retry-path leaks.

### 4. `intbig.c` `ibz_probab_prime` — trial-division precheck

`src/quaternion/ref/generic/intbig.c`

Before the Miller–Rabin loop, the candidate `n` is trial-divided by the
first 50 odd primes (3 … 233). The bulk of candidates entering MR via
`quat_represent_integer` are composite; rejecting them at `O(limbs)` per
small prime avoids paying the full `O(reps × limbs²)` Miller–Rabin cost.
Algorithm-preserving: a true prime cannot be divisible by any small prime
< itself, and any small-prime hit returns composite correctly.

### 5. `sign.c` `compute_backtracking_signature` — zero-content recovery

`src/signature/ref/lvlx/sign.c`

Return type `void → int`. When `quat_alg_make_primitive` returns a zero
`tmp` (a lattice-containment invariant violation that under `NDEBUG`
passes through silently because `quat_lattice_contains` is compiled out),
the original code performed `resp_quat->denom *= tmp = 0`, zeroing
`resp_quat->denom`. Every subsequent `quat_alg_norm` then collapsed under
`lattice_content`, `degree_full_resp / lattice_content` truncated to 0,
`random_aux_norm = 2^pow − 0` became even, and `quat_represent_integer`
fast-exited on `ibz_is_even(n_gamma)`, leaving the sign loop spinning
indefinitely until the per-iter timeout killed the process.

Fix: detect zero `tmp`, finalize locals, return 0. The caller in
`protocols_sign` treats a 0 return identically to other sub-step failures
and triggers the existing sign-retry path. Empirically this codepath was
hit in roughly 75–90 % of timeouts in our environment; it appears to be
a latent NDEBUG-masked correctness bug rather than a documented "rare
retry".

### 6. `bench.h` — per-phase wall-clock seconds (new, not in the original)

`src/common/generic/include/bench.h`

Adds a `clock_gettime(CLOCK_MONOTONIC)` pair inside the existing
`BENCH_CODE_1` / `BENCH_CODE_2` macros and emits a second printf line per
phase with median / min / max in seconds, alongside the existing
megacycles line. The paper §5.2 Table reports wall-clock seconds; the
upstream output was megacycles only, so the upstream binaries weren't
directly usable for the paper table.

### 7. `apps/benchmark.c` — `--only=<phase>` flag (new, not in the original)

`apps/benchmark.c`

Adds `--only={keypair,sign,verify,all}`. Default is `all` (preserves the
upstream behavior of timing all three phases in one process).

- `--only=keypair` times only the keypair phase; sign / verify are
  skipped.
- `--only=sign` runs keypair as untimed setup and times only sign;
  verify is skipped.
- `--only=verify` runs keypair + sign as untimed setup and times only
  verify.

Used by the bench driver so that a sign timeout cannot discard the
KeyGen / Verify samples in the same iteration.

### 8. `normeq.c` — specification-conformant `RandomIdealGivenNorm`

`src/quaternion/ref/generic/normeq.c`

Validates the documented norm/cofactor preconditions, samples the signed
coefficient intervals required by the specification, retries representation
failures, and constructs and verifies an ideal whose norm is exactly the
requested norm. Conservative bit guards reject inputs that cannot fit the
selected fixed-precision budget instead of allowing wrapped arithmetic.

### 9. `ideal.c` / `lattice.c` — bounded modular-HNF ideal construction

`src/quaternion/ref/generic/ideal.c`

Constructs `O<x,N>` in coordinates relative to `O`, reduces the generator
modulo `|N|`, and uses modular HNF bounds tied to `N`. This avoids oversized
ambient-HNF intermediates while preserving the generated ideal. Norm and HNF
invariants are checked in release builds as well as debug builds.

### 10. `intbig.c` — fixed-width safety and API correctness

`src/quaternion/ref/generic/intbig.c`

- Makes `ibz_rand_interval` inclusive, signed-endpoint aware, alias safe, and
  rejecting of invalid intervals.
- Rejects the non-representable two's-complement minimum in negation/absolute
  value paths and validates shifts, products, square roots, and primality edge
  cases.
- Uses per-call multiplication scratch rather than a shared static buffer.
- Adds optional signed-overflow checks and peak-bit tracking via
  `-DENABLE_IBZ_OVERFLOW_CHECKS=ON`.
- Enforces the paper's per-level magnitude bounds with production
  compile-time assertions and keeps zero-shift peak tracking exact.

### 11. Regression and end-to-end tests

The quaternion tests cover interval sampling, integer representation, bounded
HNF construction, exact ideal norms, and invalid inputs. The top-level
`RANDOM_IDEAL_PRECISION` test exercises protocol-sized prime and composite
norms at all three security levels and reports checked peak bit lengths. The
deterministic self-test checks KeyGen, attached and detached Sign/Verify,
tamper rejection, and in-place signing/opening.

## Build & run

```bash
cmake -B build -DSQISIGN_BUILD_TYPE=broadwell -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
ulimit -s unlimited

# Single iteration, all phases:
./build/apps/benchmark_lvl1 --iterations=1

# Per-phase timing:
./build/apps/benchmark_lvl1 --iterations=1 --only=keypair
./build/apps/benchmark_lvl1 --iterations=1 --only=sign
./build/apps/benchmark_lvl1 --iterations=1 --only=verify
```

On macOS, the self-test target requests the stack needed by the reference
implementation for levels III and V. On other platforms, `ulimit -s unlimited`
may still be required for these parameter sets.
