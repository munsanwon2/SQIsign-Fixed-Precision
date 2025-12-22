#ifndef INTBIG_H
#define INTBIG_H
#pragma once
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <rng.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C"
{
#endif

#ifndef IBZ_LIMBS
/* SQISIGN_VARIANT is defined as lvl1, lvl3, or lvl5 by the build system */
#define SQISIGN_LVL1 lvl1
#define SQISIGN_LVL3 lvl3
#define SQISIGN_LVL5 lvl5

#if SQISIGN_VARIANT == SQISIGN_LVL1
#define IBZ_LIMBS 110  /* NIST Level I */
#elif SQISIGN_VARIANT == SQISIGN_LVL3
#define IBZ_LIMBS 168  /* NIST Level III */
#elif SQISIGN_VARIANT == SQISIGN_LVL5
#define IBZ_LIMBS 222  /* NIST Level V */
#else
#define IBZ_LIMBS 222  /* Default to NIST Level V */
#endif
#endif

#ifndef IBZ_BITS
#define IBZ_BITS (IBZ_LIMBS * 64)
#endif
    typedef uint64_t ibz_t[IBZ_LIMBS];

    typedef uint64_t limb_t;
    typedef uint32_t half_limb_t;
    typedef uint64_t digit_t;

    void ibz_init(ibz_t *x);

    void ibz_finalize(ibz_t *x);


    void ibz_copy(ibz_t *target, const ibz_t *value);
    void ibz_swap(ibz_t *a, ibz_t *b);


    void ibz_neg(ibz_t *neg, const ibz_t *a);
    void ibz_abs(ibz_t *abs, const ibz_t *a);

    void ibz_add(ibz_t *sum, const ibz_t *a, const ibz_t *b);
    void ibz_sub(ibz_t *diff, const ibz_t *a, const ibz_t *b);
    void ibz_mul(ibz_t *prod, const ibz_t *a, const ibz_t *b);

    void ibz_div_2exp(ibz_t *quotient, const ibz_t *a, uint32_t exp);
    void ibz_mul_2exp(ibz_t *result, const ibz_t *a, size_t shift);


    int ibz_cmp(const ibz_t *a, const ibz_t *b);
    int ibz_is_zero(const ibz_t *x);
    int ibz_is_one(const ibz_t *x);
    int ibz_is_even(const ibz_t *x);
    int ibz_is_odd(const ibz_t *x);
    int ibz_is_negative(const ibz_t *x);

    void ibz_set(ibz_t *i, int32_t x);
    void ibz_set_u64(ibz_t *i, uint64_t x);
    int32_t ibz_get(const ibz_t *i);

    int ibz_cmp_int32(const ibz_t *x, int32_t y);
    int ibz_bitsize(const ibz_t *a);
    void ibz_to_digits(digit_t *digits, const ibz_t *a);

    /** @brief generate random value in [a, b]
     *  assumed that a >= 0 and b >= 0 and a < b
     * @returns 1 on success, 0 on failiure
     */
    int ibz_rand_interval(ibz_t *rand, const ibz_t *a, const ibz_t *b);

    /** @brief generate random value in [-m, m]
     *  assumed that m > 0 and bitlength of m < 32 bit
     * @returns 1 on success, 0 on failiure
     */
    int ibz_rand_interval_minm_m(ibz_t *rand, int32_t m);

#define ibz_to_digit_array(T, I)                                                                                       \
    do {                                                                                                               \
        memset((T), 0, sizeof(T));                                                                                     \
        ibz_to_digits((T), (I));                                                                                       \
    } while (0)
    void ibz_copy_digits(ibz_t *a, const digit_t *digits, size_t len);
#define ibz_copy_digit_array(I, T)                                                                                     \
    do {                                                                                                               \
        ibz_copy_digits((I), (T), sizeof(T) / sizeof(*(T)));                                                           \
    } while (0)
    size_t ibz_size_in_base(const ibz_t *a, int base);

    int ibz_two_adic(ibz_t *pow);

    void ibz_div(ibz_t *quotient, ibz_t *remainder, const ibz_t *a, const ibz_t *b);

    void ibz_div_floor(ibz_t *q, ibz_t *r, const ibz_t *n, const ibz_t *d);

    void ibz_mod(ibz_t *r, const ibz_t *a, const ibz_t *b);
    unsigned long ibz_mod_ui(const ibz_t *n, unsigned long d);

    void ibz_pow(ibz_t *pow, const ibz_t *x, uint32_t e);
    void ibz_pow_mod(ibz_t *pow, const ibz_t *x, const ibz_t *e, const ibz_t *m);


    void ibz_gcd(ibz_t *gcd, const ibz_t *a, const ibz_t *b);
    void ibz_gcdext(ibz_t *gcd, ibz_t *x, ibz_t *y, const ibz_t *a, const ibz_t *b);

    int ibz_invmod(ibz_t *inv, const ibz_t *a, const ibz_t *mod);

    int ibz_divides(const ibz_t *a, const ibz_t *b);

    int ibz_sqrt(ibz_t *sqrt, const ibz_t *a);

    void ibz_sqrt_floor(ibz_t *sqrt, const ibz_t *a);

    int ibz_legendre(const ibz_t *a, const ibz_t *p);

    int ibz_sqrt_mod_p(ibz_t *sqrt, const ibz_t *a, const ibz_t *p);

    extern const uint64_t ibz_const_zero[IBZ_LIMBS];
    extern const uint64_t ibz_const_one[IBZ_LIMBS];
    extern const uint64_t ibz_const_two[IBZ_LIMBS];
    extern const uint64_t ibz_const_three[IBZ_LIMBS];

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* IBZ_H */
