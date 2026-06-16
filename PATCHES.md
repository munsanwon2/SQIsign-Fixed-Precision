# Patches applied to upstream SQIsign-Fixed-Precision

This is a patched copy of the upstream **SQIsign-Fixed-Precision** baseline used
to produce the KLKL25 numbers in our bench. The patches below preserve the
original algorithm (no cryptographic change) and address build correctness,
observability, and a latent sign-loop bug. All patches are minimal.

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
token-equality `#if`. After this patch, lvl3 binaries are 168 limbs and
lvl5 binaries are 222 limbs, matching the precomp data.

### 2. `sign.c` — silence 8 `printf("--<step> done--\n")` traces

`src/signature/ref/lvlx/sign.c`

Each call to `protocols_sign` writes 8 `--<step> done--` lines to stdout.
During a sign-retry storm these can dominate stdout and drift the measured
cycle counts by roughly 10 %. Replaced each with `(void)0;`. Behavior of
the sign protocol is unchanged.

### 3. `sign.c` — `sign_retries` counter + `[SIGN_RETRIES]` stderr emit

`src/signature/ref/lvlx/sign.c`

Inside `protocols_sign`, every `continue` of the outer sign loop now
increments a `sign_retries` counter; one `fprintf(stderr, "[SIGN_RETRIES]
%d\n", sign_retries)` is emitted before `return ret`. The bench driver
reads this so a single sign that took many retries can be distinguished
from a sign that took one fast attempt. Observability only.

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

`ulimit -s unlimited` is required (especially at lvl3 / lvl5). Sign
retries are visible in stderr as `[SIGN_RETRIES] <N>`.
