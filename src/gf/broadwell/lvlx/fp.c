// SPDX-License-Identifier: Apache-2.0
// Adapted from the SQIsign nist-v3 saturated field support (6d017708db403bf83977fa70770fc4f7f9e9ff21).
// The field elements use R = 2^(64 * NWORDS_FIELD) and are below 2^FP_FIELD_BITS.
#include <encoded_sizes.h>
#include <fp.h>
#include <string.h>

#define N NWORDS_FIELD
#define TOPW (FP_FIELD_BITS - 64 * (N - 1)) // bits used in the top limb

__extension__ typedef unsigned __int128 fpg_u128_t;

const digit_t fp_modulus[N] = FP_MODULUS_INIT;
const digit_t fp_twice_modulus[N] = FP_TWICE_MODULUS_INIT;
static const fp_t R2 = { FP_MONT_R2 };
static const fp_t TWO_INV = { FP_MONT_TWO_INV };
static const fp_t THREE_INV = { FP_MONT_THREE_INV };
static const fp_t ONE_RAW = { { 1 } };
const fp_t ZERO = { { 0 } };
const fp_t ONE = { FP_MONT_ONE };

static uint64_t
mask64(uint32_t ctl)
{
    return (uint64_t)(int64_t)(int32_t)ctl;
}

static uint32_t
limbs_iszero(const uint64_t *v, const uint64_t *xorwith)
{
    uint64_t r = 0;
    for (int i = 0; i < N; i++)
        r |= v[i] ^ (xorwith ? xorwith[i] : 0);
    return (uint32_t)((r | (0 - r)) >> 63) - 1;
}

// w (N+1 limbs, < 2^(FP_TORSION_POWER+64)) -> d in [0, 2^FP_FIELD_BITS), value preserved mod p
static void
partial_reduce(fp_t *d, const uint64_t *w)
{
    const int ti = FP_TORSION_POWER / 64, ts = FP_TORSION_POWER % 64;
    uint64_t h = w[ti] >> ts;
    if (ts)
        h |= w[ti + 1] << (64 - ts);
    uint64_t q = h / FP_COFACTOR, r = h - q * FP_COFACTOR;

    uint64_t out[NWORDS_FIELD];
    for (int i = 0; i < N; i++)
        out[i] = (i < ti) ? w[i] : 0;
    if (ts)
        out[ti] = w[ti] & ((UINT64_C(1) << ts) - 1);

    // out += q + (r << FP_TORSION_POWER)
    fpg_u128_t cc = (fpg_u128_t)out[0] + q;
    out[0] = (uint64_t)cc;
    cc >>= 64;
    for (int i = 1; i < N; i++) {
        cc += out[i];
        out[i] = (uint64_t)cc;
        cc >>= 64;
    }
    fpg_u128_t rs = (fpg_u128_t)r << ts;
    cc = (fpg_u128_t)out[ti] + (uint64_t)rs;
    out[ti] = (uint64_t)cc;
    cc >>= 64;
    cc += (uint64_t)(rs >> 64);
    for (int i = ti + 1; i < N; i++) {
        cc += out[i];
        out[i] = (uint64_t)cc;
        cc >>= 64;
    }

    // one conditional subtraction of p (result <= p + q < 2^BITS + 2^64)
    uint64_t k = 0 - (out[N - 1] >> TOPW);
    unsigned char b = 0;
    for (int i = 0; i < N; i++) {
        fpg_u128_t t = (fpg_u128_t)out[i] - (fp_modulus[i] & k) - b;
        d->v[i] = (uint64_t)t;
        b = (t >> 64) & 1;
    }
}

void
fp_copy(fp_t *d, const fp_t *a)
{
    memcpy(d, a, sizeof(fp_t));
}

void
fp_set_zero(fp_t *a)
{
    memcpy(a, &ZERO, sizeof(fp_t));
}

void
fp_set_one(fp_t *a)
{
    memcpy(a, &ONE, sizeof(fp_t));
}

void
fp_neg(fp_t *d, const fp_t *a)
{
    fp_sub(d, &ZERO, a);
}

void
fp_half(fp_t *d, const fp_t *a)
{
    fp_mul(d, a, &TWO_INV);
}

void
fp_div3(fp_t *d, const fp_t *a)
{
    fp_mul(d, a, &THREE_INV);
}

void
fp_set_small(fp_t *d, digit_t x)
{
    fp_t t = { { 0 } };
    t.v[0] = x;
    fp_mul(d, &t, &R2);
}

void
fp_mul_small(fp_t *d, const fp_t *a, uint32_t x)
{
    uint64_t w[NWORDS_FIELD + 1];
    fpg_u128_t cc = 0;
    for (int i = 0; i < N; i++) {
        cc += (fpg_u128_t)a->v[i] * x;
        w[i] = (uint64_t)cc;
        cc >>= 64;
    }
    w[N] = (uint64_t)cc;
    partial_reduce(d, w);
}

void
fp_select(fp_t *d, const fp_t *a0, const fp_t *a1, uint32_t ctl)
{
    uint64_t m = mask64(ctl);
    for (int i = 0; i < N; i++)
        d->v[i] = a0->v[i] ^ (m & (a0->v[i] ^ a1->v[i]));
}

void
fp_cswap(fp_t *a, fp_t *b, uint32_t ctl)
{
    uint64_t m = mask64(ctl);
    for (int i = 0; i < N; i++) {
        uint64_t t = m & (a->v[i] ^ b->v[i]);
        a->v[i] ^= t;
        b->v[i] ^= t;
    }
}

uint32_t
fp_is_zero(const fp_t *a)
{
    return limbs_iszero(a->v, NULL) | limbs_iszero(a->v, fp_modulus);
}

uint32_t
fp_is_equal(const fp_t *a, const fp_t *b)
{
    fp_t t;
    fp_sub(&t, a, b);
    return fp_is_zero(&t);
}

// canonical value of a into m (in [0, p))
static void
canonical(uint64_t *m, const fp_t *a)
{
    fp_t t;
    fp_mul(&t, a, &ONE_RAW); // leave Montgomery; result < 2^BITS < 2p
    unsigned char b = 0;
    uint64_t s[NWORDS_FIELD];
    for (int i = 0; i < N; i++) {
        fpg_u128_t u = (fpg_u128_t)t.v[i] - fp_modulus[i] - b;
        s[i] = (uint64_t)u;
        b = (u >> 64) & 1;
    }
    uint64_t keep = 0 - (uint64_t)b; // borrow -> t < p -> keep t
    for (int i = 0; i < N; i++)
        m[i] = (t.v[i] & keep) | (s[i] & ~keep);
}

_Static_assert(FP_ENCODED_BYTES == 8 * N && N * 64 - FP_FIELD_BITS >= 4,
               "unsupported saturated field layout");

void
fp_encode(void *dst, const fp_t *a)
{
    uint64_t m[N];
    canonical(m, a);
    uint8_t *out = dst;
    for (size_t i = 0; i < FP_ENCODED_BYTES; i++)
        out[i] = (uint8_t)(m[i / 8] >> (8 * (i % 8)));
}

uint32_t
fp_decode(fp_t *d, const void *src)
{
    const uint8_t *in = src;
    fp_t v = { { 0 } };
    for (size_t i = 0; i < FP_ENCODED_BYTES; i++)
        v.v[i / 8] |= (uint64_t)in[i] << (8 * (i % 8));

    unsigned char borrow = 0;
    for (size_t i = 0; i < N; i++) {
        fpg_u128_t u = (fpg_u128_t)v.v[i] - fp_modulus[i] - borrow;
        borrow = (u >> 64) & 1;
    }
    uint64_t mask = (uint64_t)0 - borrow;
    for (size_t i = 0; i < N; i++)
        v.v[i] &= mask;
    fp_mul(d, &v, &R2);
    return (uint32_t)mask;
}

void
fp_decode_reduce(fp_t *d, const void *src, size_t len)
{
    const uint8_t *s = src;
    const size_t B = 8 * (size_t)N;
    fp_t acc = { { 0 } };
    if (len == 0) {
        *d = acc;
        return;
    }
    size_t nb = (len + B - 1) / B;
    for (size_t blk = nb; blk-- > 0;) {
        size_t off = blk * B;
        size_t n = (off + B <= len) ? B : len - off;
        uint64_t w[NWORDS_FIELD + 1];
        memset(w, 0, sizeof w);
        for (size_t i = 0; i < n; i++)
            w[i / 8] |= (uint64_t)s[off + i] << (8 * (i % 8));
        fp_t bred;
        partial_reduce(&bred, w);
        if (blk == nb - 1) {
            acc = bred;
        } else {
            fp_t sh;
            fp_mul(&sh, &acc, &R2); // acc * 2^(64N) mod p
            fp_add(&acc, &sh, &bred);
        }
    }
    fp_mul(d, &acc, &R2); // to Montgomery
}

// d = a^(2^k - 1), k >= 1
static void
pow_ones(fp_t *d, const fp_t *a, int k)
{
    fp_t x = *a;
    int e = 1;
    for (int i = 62 - __builtin_clzll((unsigned long long)k); i >= 0; i--) {
        fp_t y = x;
        for (int j = 0; j < e; j++)
            fp_sqr(&y, &y);
        fp_mul(&x, &y, &x);
        e <<= 1;
        if ((k >> i) & 1) {
            fp_sqr(&x, &x);
            fp_mul(&x, &x, a);
            e += 1;
        }
    }
    *d = x;
}

// out <- a^((p-3)/4), with (p-3)/4 = (C-1)*2^(T-2) + (2^(T-2) - 1)
void
fp_exp3div4(fp_t *out, const fp_t *a)
{
    fp_t base = *a, M;
    pow_ones(&M, &base, FP_TORSION_POWER - 2);
#if FP_COFACTOR == 1
    *out = M;
#else
    fp_t Ma, z;
    fp_mul(&Ma, &M, &base); // a^(2^(T-2))
    z = Ma;
    for (int i = 62 - __builtin_clzll((unsigned long long)(FP_COFACTOR - 1)); i >= 0; i--) {
        fp_sqr(&z, &z);
        if (((FP_COFACTOR - 1) >> i) & 1)
            fp_mul(&z, &z, &Ma);
    }
    fp_mul(out, &z, &M);
#endif
}

void
fp_inv(fp_t *a)
{
    fp_t a0 = *a, t = *a;
    fp_exp3div4(&t, &t);
    fp_sqr(&t, &t);
    fp_sqr(&t, &t);
    fp_mul(a, &t, &a0); // a^(p-2)
}

void
fp_sqrt(fp_t *a)
{
    fp_t a0 = *a, t = *a;
    fp_exp3div4(&t, &t);
    fp_mul(a, &t, &a0); // a^((p+1)/4)
}

uint32_t
fp_is_square(const fp_t *a)
{
    fp_t t = *a;
    fp_exp3div4(&t, &t);
    fp_sqr(&t, &t);
    fp_mul(&t, &t, a); // a^((p-1)/2)
    return fp_is_equal(&t, &ONE) | fp_is_zero(a);
}
