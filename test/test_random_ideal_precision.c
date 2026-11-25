// SPDX-License-Identifier: Apache-2.0

#include <stdio.h>

#include <encoded_sizes.h>
#include <quaternion.h>
#include <quaternion_constants.h>
#include <quaternion_data.h>
#include <rng.h>
#include <torsion_constants.h>

#include "internal.h"

#ifndef EXPECTED_IBZ_LIMBS
#error "EXPECTED_IBZ_LIMBS must be supplied by CMake"
#endif

_Static_assert(IBZ_LIMBS == EXPECTED_IBZ_LIMBS, "wrong fixed-precision budget for this variant");
_Static_assert(IBZ_BITS - 1 >= PAPER_IBZ_BITS && IBZ_BITS - 1 >= LINE13_IBZ_BITS,
               "signed fixed-precision capacity is below a Section 3 expression");

static int
sample_and_check(const char *label,
                 const ibz_t *norm,
                 int is_prime,
                 const ibz_t *cofactor,
                 const quat_represent_integer_params_t *params)
{
    quat_left_ideal_t ideal;
    int ok;

    quat_left_ideal_init(&ideal);

#ifdef IBZ_ENABLE_OVERFLOW_CHECKS
    ibz_overflow_checks_reset_peak();
#endif

    printf("%s: sampling...\n", label);
    fflush(stdout);
    ok = quat_sampling_random_ideal_O0_given_norm(&ideal, norm, is_prime, params, cofactor);
    if (!ok || ibz_cmp(&ideal.norm, norm) != 0 || ideal.parent_order != &params->order->order ||
        !ibz_mat_4x4_is_hnf(&ideal.lattice.basis) ||
        !quat_lattice_inclusion(&ideal.lattice, ideal.parent_order)) {
        fprintf(stderr, "%s sampling failed or returned an invalid ideal\n", label);
        quat_left_ideal_finalize(&ideal);
        return 0;
    }

    // Recompute the norm from the lattice, rather than trusting the cached field.
    printf("%s: recomputing lattice norm...\n", label);
    fflush(stdout);
    quat_lideal_norm(&ideal);
    if (ibz_cmp(&ideal.norm, norm) != 0) {
        fprintf(stderr, "%s recomputed norm mismatch\n", label);
        quat_left_ideal_finalize(&ideal);
        return 0;
    }

    const int n_bits = ibz_bitsize(norm);
    const int p_bits = ibz_bitsize(&params->algebra->p);
    const int analytic_hnf_bound = 8 * n_bits + 4 * p_bits + 13;
    if (analytic_hnf_bound > IBZ_BITS - 1) {
        fprintf(stderr, "%s analytic HNF bound exceeds the signed precision budget\n", label);
        quat_left_ideal_finalize(&ideal);
        return 0;
    }

#ifdef IBZ_ENABLE_OVERFLOW_CHECKS
    const int peak = ibz_overflow_checks_get_peak();
    if (peak > IBZ_BITS - 1) {
        fprintf(stderr, "%s measured peak exceeds signed capacity\n", label);
        quat_left_ideal_finalize(&ideal);
        return 0;
    }
    printf("%s: N=%d bits, p=%d bits, analytic-bound=%d, checked-peak=%d, capacity=%d\n",
           label,
           n_bits,
           p_bits,
           analytic_hnf_bound,
           peak,
           IBZ_BITS - 1);
#else
    printf("%s: N=%d bits, p=%d bits, analytic-bound=%d, capacity=%d\n",
           label,
           n_bits,
           p_bits,
           analytic_hnf_bound,
           IBZ_BITS - 1);
#endif

    quat_left_ideal_finalize(&ideal);
    return 1;
}

static void
maximum_auxiliary_composite(ibz_t *norm, const ibz_t *cofactor)
{
    ibz_t factor, gcd;
    ibz_init(&factor);
    ibz_init(&gcd);

    // A large, explicitly composite odd representative within the protocol's
    // response-length budget: N = 3*(2^(r-2)-1), adjusted only for coprimality.
    ibz_pow(&factor, &ibz_const_two, SQIsign_response_length - 2);
    ibz_sub(&factor, &factor, &ibz_const_one);
    do {
        ibz_mul(norm, &factor, &ibz_const_three);
        ibz_gcd(&gcd, norm, cofactor);
        if (!ibz_is_one(&gcd)) {
            ibz_sub(&factor, &factor, &ibz_const_two);
        }
    } while (!ibz_is_one(&gcd));

    ibz_finalize(&gcd);
    ibz_finalize(&factor);
}

int
main(void)
{
    unsigned char seed[48] = { 0 };
    ibz_t auxiliary_norm;
    int result = 0;

    seed[0] = (unsigned char)EXPECTED_IBZ_LIMBS;
    seed[1] = EXPECTED_IBZ_LIMBS >> 8;
    seed[2] = 0x52;
    seed[3] = 0x49;
    randombytes_init(seed, NULL, 256);

    const quat_represent_integer_params_t params = {
        .algebra = &QUATALG_PINFTY,
        .order = &EXTREMAL_ORDERS[0],
        .primality_test_iterations = QUAT_primality_num_iter,
    };

    printf("RandomIdealGivenNorm fixed-precision test: limbs=%d, signed-capacity=%d bits\n",
           IBZ_LIMBS,
           IBZ_BITS - 1);

    if (ibz_bitsize(&SEC_DEGREE) != EXPECTED_MIX_BITS || ibz_cmp(&SEC_DEGREE, &COM_DEGREE) != 0) {
        fprintf(stderr, "invalid mixing degree for this variant\n");
        return 1;
    }
    result |= !sample_and_check("SEC_DEGREE", &SEC_DEGREE, 1, NULL, &params);
    result |= !sample_and_check("COM_DEGREE", &COM_DEGREE, 1, NULL, &params);

    ibz_init(&auxiliary_norm);
    maximum_auxiliary_composite(&auxiliary_norm, &QUAT_prime_cofactor);
    result |= !sample_and_check(
        "maximum auxiliary composite", &auxiliary_norm, 0, &QUAT_prime_cofactor, &params);
    ibz_finalize(&auxiliary_norm);

    return result;
}
