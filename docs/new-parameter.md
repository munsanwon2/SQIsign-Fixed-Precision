# Round-2 algorithms with Round-3 prime/security parameters

This branch starts at `d0cb037ee6a9f68f55cb6f55b4e3746c79550330`. It uses the
three main primes, explicit security levels and mixing degrees from SQIsign v3.0
with the Round-2 quaternion, Clapotis and signature algorithms. Field kernels
come from the Apache-2.0 upstream `nist-v3` tag at
`6d017708db403bf83977fa70770fc4f7f9e9ff21`.
[Round-3 specification](https://sqisign.org/spec/sqisign-20260901.pdf)
[Round-2 specification](https://sqisign.org/spec/sqisign-20250707.pdf)
[Upstream source](https://github.com/SQISign/the-sqisign/tree/nist-v3)

| Level | Prime | λ | f | `D_mix` bits | Response exponent | Challenge space / scalar bits | Hash iterations |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| I | 3 × 2^324 − 1 | 128 | 324 | 582 | 163 | 161 / 128 | 1 |
| III | 27 × 2^500 − 1 | 192 | 500 | 889 | 253 | 247 / 192 | 1 |
| V | 17 × 2^664 − 1 | 256 | 664 | 1181 | 335 | 329 / 256 | 1 |

`SEC_DEGREE` and `COM_DEGREE` both equal the first prime after
`p · 2^(2λ)`. The challenge isogeny continues to use the full `2^f` torsion.
The Round-2 response exponent is `ceil(bit_length(p)/2)`. Verification uses
up to that exponent plus two extra HD torsion levels, which fits within `f`
for every variant. The SHAKE construction still absorbs the two encoded
j-invariants and message, then hashes the intermediate digest again before
clamping the output. With `HASH_ITERATIONS=1` it performs two SHAKE operations
and zero loop iterations.

| Level | Field / scalar storage bytes | Secret-ideal entry bytes | Response matrix entry bytes | Each response counter, bytes | PK / SK / signature bytes |
| --- | ---: | ---: | ---: | ---: | ---: |
| I | 48 / 48 | 48 | 21 | 1 | 97 / 501 / 200 |
| III | 64 / 64 | 66 | 32 | 1 | 129 / 711 / 284 |
| V | 88 / 88 | 88 | 43 | 2 | 177 / 953 / 386 |

Field encodings are little-endian and padded to a multiple of eight bytes.
Field decoders reject integers at or above `p` and nonzero padding. The
reference backend uses unsaturated Montgomery limbs (64-bit word radices
55/57/61 and 32-bit word radices 28/29/28). Broadwell uses saturated 64-bit
Montgomery limbs. Both export the same field values and wire format. Scalar
words follow the selected radix and are converted explicitly to and from the
64-bit-limb quaternion integer format.

Counters use little-endian encoding and Level V takes two bytes per counter
to cover the 335-bit response. Verification checks the combined length
before its matrix-bound shift. Level III's five secret-ideal entries use 66
bytes each. The equivalent-prime norm bound is below `2^19 p`, and 512 bits
would only cover part of this range. For a prime norm `N`, at least one HNF
column or a sum of two columns has norm divisible by `N` exactly once.
Otherwise polarization would make all norms in the ideal divisible by
`N²`. Since `N O0 ⊂ I`, the HNF row pivots are at most `2N` in the
integral numerator basis. Four times the norm bounds the generator
coordinates, including either sign. The `SK` encoding retains the public key
prefix and the four full-torsion scalar entries.

## Section 3 precision

| Level | KeyGen expression | Paper magnitude bits | Line-13 expression bits | Larger verified expression | 64-bit limbs | Signed capacity |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| I | 2,043 | 9,555 | 9,565 | 9,565 | 150 | 9,599 |
| III | 3,118 | 14,572 | 14,569 | 14,572 | 228 | 14,591 |
| V | 4,098 | 19,145 | 19,145 | 19,145 | 300 | 19,199 |

The task gives the expressions and Table 3 values for the manuscript with
SHA-256 `94d946d798f35798b84aa5de6b69adf2024e0487720f292523d0f6ee250c540c`.
The task document was available and the published ePrint PDF has a different
hash and earlier budgets. The parameter audit recomputes the expressions
from the supplied task with 256-bit real interval arithmetic and checks the
ceilings at both endpoints. [Public ePrint record](https://eprint.iacr.org/2026/1031)

The ratio of the preceding line-13 expression to the paper's final Sign
expression is `2^16 (2^f / p)^4`. It is greater than one at Level I. The
production assertions check the original and preceding expressions against
the positive two's-complement capacity. The manifest records both values
and leaves `validated_required_bits` empty while the full code-level bound
remains open.

### Intersection code and bound coverage

Signing first intersects the secret and challenge ideals and then intersects
that result with the row-sign conjugate of the commitment ideal. Integral
ideals in `O0` have numerator denominator at most two. For an ideal of norm
`N` with denominator `d`, the numerator matrix determinant has magnitude
`D=d^4 N²/4`, using `[O0:I]=N²`. The first intersection has norm
`2^f N_sk` because the two indices are coprime. The last operand in the
second intersection retains the absolute row-pivot bounds after conjugation.

`ibz_mat_4x4_is_hnf` and the HNF column reductions use an upper-triangular
column basis. Every entry in row `i` has magnitude at most the pivot `h_i`.
The 2-by-2 minors, including their partial products and differences, have
magnitude at most twice the product of the two row pivots. Each adjugate
partial sum is bounded by six times the product of the three relevant row
pivots. The determinant loop adds up to six products of two minors, giving a
bound of `24D` for its products and partial sums. This argument uses the
actual minor expansion and applies to the signed-row variant.

In `quat_lattice_add`, the cross-multiplied dual for input `i` has row bound
`6 d_i D_1 D_2 / h_i`. Its determinant-magnitude bound, including the
factors from the minors and partial sums, is
`24 · 6^4 · d_i^4 · (D_1 D_2)^4 / D_i`. The HNF modulus before denominator
reduction is exactly
`(D_1 D_2)^3 · (d_1 d_2)^4 · gcd(N_1,N_2)^2 / 4`. The first pair has gcd
one. For the second pair the bound allows the commitment and secret norms to
coincide. Using `N_sk,N_com < 2^19 p` gives the following conservative ceilings.

| Level | Scaled determinant intermediates, first / second | Add HNF modulus, first / second | Final HNF modulus, first / second |
| --- | ---: | ---: | ---: |
| I | 4,737 / 7,455 | 4,032 / 6,792 | 1,340 / 2,030 |
| III | 7,225 / 11,369 | 6,162 / 10,354 | 2,050 / 3,098 |
| V | 9,521 / 14,977 | 8,130 / 13,634 | 2,706 / 4,082 |

The reduced denominator of the sum of duals divides `lcm(N_1,N_2)`
because each original ideal contains `N_i O0`. The dual of that sum is the
intersection, which lies in `O0`. The existing `quat_lattice_hnf` reduces
the denominator before the final modular HNF. The resulting determinant,
and hence the last HNF modulus, is at most `4 (N_1 N_2)^2`. The lower-triangular
matrix at that point inherits row bounds from the adjugate of the dual-sum
HNF. Applying the same minor expansion bounds its determinant products and
partial sums by `24 · 6^4` times the determinant magnitude.

The remaining gap is a capacity bound for all exact Bezout products and
column reductions *inside* the HNF of the cross-multiplied dual sum.
Squaring the conservative second-call modulus exceeds the storage budget
and a tighter pivot-dependent argument would be needed. The overflow
checks exercise the actual path, including its HNF temporaries. Their
measured peaks supply execution evidence. The branch keeps the existing
intersection and lattice reduction algorithms.

The generator checks maximality and odd integral connecting-ideal norms.
The largest alternate-order `q` values are 113, 181 and 61. The audit checks
`norm² < q p` for each generated connecting ideal and recomputes the signed
initializers, including zero, negative values and 64-bit boundaries. Sage
checks `p` and the next-prime calculation with proofs enabled and uses
full-order pairings to check the torsion and the endomorphism action.

## Regeneration and tests

Use SageMath 10.10, the `Theta-SageMath` directory from
`ThetaIsogenies/two-isogenies` at
`0ea5020effabde4bdc5e6f425ba5b94209f0e58b` and `Jonathke/deuring-2D` at
`799d69d0c784a57446566e15be4207b936bcc469`. The pinned Deuring revision
provides the `SuitableIdeals` interface used by the generator.
The endomorphism generators can take substantially longer for Level V.

```sh
export PYTHONPATH="/path/to/deuring-2D:/path/to/two-isogenies/Theta-SageMath"
cmake -S . -B build-precomp
cmake --build build-precomp --target precomp -j1
python3 scripts/precomp/make_manifest.py
python3 scripts/precomp/import_field_kernels.py /path/to/the-sqisign
```

The repository includes the generated C and headers. A C build needs no
SageMath. Use an expanded stack for signing and full protocol tests. The
largest FindUV arrays are roughly 270 MB in Level V after removing the two
unused arrays. The threaded test reserves 512 MiB per thread.

```sh
ulimit -s unlimited
cmake -S . -B build-debug -DSQISIGN_BUILD_TYPE=ref -DCMAKE_BUILD_TYPE=Debug \
  -DENABLE_GMP=OFF -DENABLE_IBZ_OVERFLOW_CHECKS=ON
cmake --build build-debug -j2
ctest --test-dir build-debug --output-on-failure -R 'SELFTEST|RANDOM_IDEAL_PRECISION|PARAMETER_ENCODING'

cmake -S . -B build-release -DSQISIGN_BUILD_TYPE=ref -DCMAKE_BUILD_TYPE=Release -DENABLE_GMP=OFF
cmake --build build-release -j2
scripts/gen_kat_files.sh build-release --count=2
ctest --test-dir build-release --output-on-failure -R 'SELFTEST|KAT'
```

The branch fixtures use the first two cases of the repository's
`PQCgenKAT_sign` generator. Its 48-byte initial entropy is `00 01 … 2f`.
The test RNG uses NIST AES-256 CTR-DRBG and reseeds from each fixture's
48-byte seed. The signed response is verified independently by the API,
and the reader regenerates the keys and signatures byte for byte.
[Parameter manifest](new-parameter-manifest.json)

## Verification results

The test host was x86-64 with BMI2, ADX, AES and AVX2, GCC 14.2.0 and
CMake 4.4.5. The KAT generator source has Git blob `8413e77b08cd32a91b956bf4b937084015b5900c`.

- The overflow-check reference Debug suite passed 24 tests for the field,
  curve and pairing code, quaternion arithmetic, random ideals and wire
  encodings. The quaternion tests also passed after the aliased-division
  regression check was added.
- The overflow-check reference self-tests passed at all three levels with
  `--seed=1 --msglen=32`. They cover attached and detached verification,
  in-place signing and opening, and rejection of changed messages and
  signatures. The signing and verification paths exercise the HD chains.
- `ctest -R 'KAT|PARAMETER_ENCODING'` reproduced both KAT vectors per
  level in the 64-bit reference, 32-bit `opt` and Broadwell Release builds.
  The 32-bit build sets `-DGF_RADIX=32 -DSQISIGN_BUILD_TYPE=opt` on the
  same x86-64 host. The Broadwell build sets `-DSQISIGN_BUILD_TYPE=broadwell`.
  All three builds set `-DENABLE_GMP=OFF` and their KAT binaries have no
  GMP dependency.
- A Python integer oracle matched 37 field-vector pairs on each of the nine
  level/backend combinations. It checked padded encodings, addition,
  subtraction, multiplication, squaring, inversion, square roots and
  reduction of double-width inputs.
- A Python integer oracle matched all 1,208 signed division cases up to
  19,199 magnitude bits, including aliased input and output buffers. The
  division test library ran with UBSan and fixed-precision overflow checks.
- The verification-only ASan/UBSan build passed the three field tests and
  the three wire-encoding tests and verified both KAT signatures at each
  level. The malformed cases include out-of-range response counters,
  shortened inputs and noncanonical field values.
- A C audit checked the full torsion order, the quaternion relations on the
  torsion matrices, the Frobenius action and the connecting-ideal Weil
  pairing relation for all 22 starting curves. Separate Sage runs reproduced
  all nine Level I precomputation steps and the Level III and V outputs
  from eight steps in isolated directories. Level III's remaining
  endomorphism tables also matched a separate calculation using another
  discrete-log procedure. A second full Level V endomorphism run remains
  open.

| Level | Observed self-test peak, magnitude bits | Signed capacity |
| --- | ---: | ---: |
| I | 6,756 | 9,599 |
| III | 10,291 | 14,591 |
| V | 14,226 | 19,199 |

The isolated ID-to-isogeny stress suite at Level I reached the 1,500-second
test timeout. The full stress suites for Levels III and V and the ARM
Cortex-M4 hardware path were outside the completed test runs.
