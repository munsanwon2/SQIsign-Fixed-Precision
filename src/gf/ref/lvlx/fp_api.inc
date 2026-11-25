/* SPDX-License-Identifier: Apache-2.0 */

#include <encoded_sizes.h>
#include <fp.h>

_Static_assert(Wordlength == RADIX && Radix == FP_LIMB_RADIX && Nlimbs == NWORDS_FIELD &&
                   Nbits == FP_FIELD_BITS && Nbytes == (Nbits + 7) / 8,
               "field kernel and precomputation use different parameters");
_Static_assert(Nbytes <= FP_ENCODED_BYTES, "field encoding is too short");

const digit_t ZERO[NWORDS_FIELD] = { 0 };
const digit_t ONE[NWORDS_FIELD] = FP_MONT_ONE;
static const fp_t TWO_INV = FP_MONT_TWO_INV;
static const fp_t THREE_INV = FP_MONT_THREE_INV;

void
fp_set_small(fp_t *x, const digit_t val)
{
    const spint mask = ((spint)1 << Radix) - 1;
    const spint limbs[Nlimbs] = { val & mask, val >> Radix };
    nres(limbs, *x);
}

void
fp_mul_small(fp_t *x, const fp_t *a, const uint32_t val)
{
    fp_t b;
    fp_set_small(&b, val);
    fp_mul(x, a, &b);
}

void
fp_set_zero(fp_t *x)
{
    modzer(*x);
}

void
fp_set_one(fp_t *x)
{
    modone(*x);
}

uint32_t
fp_is_equal(const fp_t *a, const fp_t *b)
{
    return -(uint32_t)modcmp(*a, *b);
}

uint32_t
fp_is_zero(const fp_t *a)
{
    return -(uint32_t)modis0(*a);
}

void
fp_copy(fp_t *out, const fp_t *a)
{
    modcpy(*a, *out);
}

void
fp_cswap(fp_t *a, fp_t *b, uint32_t ctl)
{
    digit_t mask = (digit_t)(0 - (digit_t)(ctl & 1));
    for (size_t i = 0; i < Nlimbs; i++) {
        digit_t tmp = mask & ((*a)[i] ^ (*b)[i]);
        (*a)[i] ^= tmp;
        (*b)[i] ^= tmp;
    }
}

void
fp_add(fp_t *out, const fp_t *a, const fp_t *b)
{
    modadd(*a, *b, *out);
}

void
fp_sub(fp_t *out, const fp_t *a, const fp_t *b)
{
    modsub(*a, *b, *out);
}

void
fp_neg(fp_t *out, const fp_t *a)
{
    modneg(*a, *out);
}

void
fp_sqr(fp_t *out, const fp_t *a)
{
    modsqr(*a, *out);
}

void
fp_mul(fp_t *out, const fp_t *a, const fp_t *b)
{
    modmul(*a, *b, *out);
}

void
fp_inv(fp_t *x)
{
    modinv(*x, NULL, *x);
}

uint32_t
fp_is_square(const fp_t *a)
{
    return -(uint32_t)modqr(NULL, *a);
}

void
fp_sqrt(fp_t *a)
{
    modsqrt(*a, NULL, *a);
}

void
fp_half(fp_t *out, const fp_t *a)
{
    modmul(TWO_INV, *a, *out);
}

void
fp_exp3div4(fp_t *out, const fp_t *a)
{
    modpro(*a, *out);
}

void
fp_div3(fp_t *out, const fp_t *a)
{
    modmul(THREE_INV, *a, *out);
}

void
fp_encode(void *dst, const fp_t *a)
{
    unsigned char *out = dst;
    spint normal[Nlimbs];
    redc(*a, normal);
    for (size_t i = 0; i < Nbytes; i++) {
        out[i] = (unsigned char)normal[0];
        (void)modshr(8, normal);
    }
    memset(out + Nbytes, 0, FP_ENCODED_BYTES - Nbytes);
}

uint32_t
fp_decode(fp_t *d, const void *src)
{
    const unsigned char *in = src;
    const unsigned int high_mask = (1U << (Nbits - 8 * (Nbytes - 1))) - 1;
    uint32_t excess = in[Nbytes - 1] & ~high_mask;
    for (size_t i = Nbytes; i < FP_ENCODED_BYTES; i++)
        excess |= in[i];

    spint normal[Nlimbs] = { in[Nbytes - 1] & high_mask };
    for (size_t i = Nbytes - 1; i-- > 0;) {
        modshl(8, normal);
        normal[0] |= in[i];
    }
    uint32_t canonical = (uint32_t)modfsb(normal);
    uint32_t ok = ((excess - 1) >> 31) & canonical;
    spint mask = (spint)0 - ok;
    for (size_t i = 0; i < Nlimbs; i++)
        normal[i] &= mask;
    nres(normal, *d);
    return 0U - ok;
}

void
fp_decode_reduce(fp_t *d, const void *src, size_t len)
{
    const unsigned char *in = src;
    fp_set_zero(d);
    for (size_t i = len; i-- > 0;) {
        fp_t byte;
        for (int bit = 0; bit < 8; bit++)
            fp_add(d, d, d);
        modint(in[i], byte);
        fp_add(d, d, &byte);
    }
}
